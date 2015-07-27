classdef modelFitter < handle
    %MODELFITTER GUI class for maximum likelihood fitting of single channel
    %mechanism
        
    properties
        %GUI Object Handles
        Figure
        OpenTimeAxis
        ShutTimeAxis
        RunBtn
        SetIdlDataBtn
        IdlDataFileLbl
        OpenResolution
        ShutResolution
        ImposeResBtn
        BurstsChkBox
        TCrit
        ModelSelectBtn
        ModelEditBtn
        
        %GUI Data
        IdlFilename
        Dwells
        States
        ResolvedDwells
        ResolvedStates
        Bursts
        BurstStates
        ModelFilename
        Model = gatingmodel;
    end
    
    methods
        
        function app = timeCourseFitter
            scrsz = get(0,'ScreenSize');
            pos = [50+0.05*scrsz(1), 75+0.05*scrsz(2), 0.85*scrsz(3), 0.75*scrsz(4)];
            app.Figure = figure('MenuBar','none','NumberTitle','off',...
                'Name','Single Channel Model Fitting','Toolbar','figure',...
                'Position',pos,'CloseRequestFcn',@app.closeApp,...
                'Visible','on');
            app.OpenTimeAxis = axes('Parent',app.Figure,...
                'Position',[.05 .1 .9 .75]);
            x=0.1*pos(3);
            y=0.92*pos(4);
            wid=100;
            hgt=30;
            dx = 10;
            app.OpenBtn = uicontrol(app.Figure,'Style','push',...
                'String','Open Data','Position',[x y wid hgt],...
                'Callback',@app.openData);
            app.FitBtn = uicontrol(app.Figure,'Style','push',...
                'String','Time Course Fitting','Callback',@app.fitData,...
                'Position',[x+wid+dx y wid hgt]);
            app.FilterBtn = uicontrol(app.Figure,'Style','toggle',...
                'String','Filter','Position',[x+2*(wid+dx) y wid hgt],...
                'Value',0,'Callback',@app.filterData);
            app.FilterValTxt = uicontrol(app.Figure,'Style','edit',...
                'Position',[x+3*(wid+dx) y wid hgt],'Callback',@app.filterVal,...
                'ButtonDownFcn',@app.filterValBtnDwn,'BackgroundColor','w');
            filterLabel = uicontrol(app.Figure,'Style','text','String','kHz',...
                'Position',[x+4*(wid+dx) y+round(hgt/4) 30 round(hgt/2)],...
                'HorizontalAlign','left');
            app.AmpHistBtn = uicontrol(app.Figure,'Style','push',...
                'String','Amplitude Histogram','Callback',@app.ampHist,...
                'Position',[x+4*(wid+dx)+30+dx y wid hgt]);
            app.BaseSubBtn = uicontrol(app.Figure,'Style','push',...
                'String','Subtract Baseline','Callback',@app.baselineCorrect,...
                'Position',get(app.AmpHistBtn,'Position')+[dx+wid 0 0 0]);
            app.FileMenu = uimenu(app.Figure,'Label','File');
                app.OpenMenu = uimenu(app.FileMenu,'Label','Open ABF',...
                    'Callback',@app.openData);
                uimenu(app.FileMenu,'Label','Export To Workspace',...
                    'callback',@app.exportToWksp);
                uimenu(app.FileMenu,'Label','Import From Workspace',...
                    'callback',@app.importFromWksp);
                app.SaveIdlMenu = uimenu(app.FileMenu,'Label','Save Idealization',...
                    'Callback',@app.saveIdl);
                uimenu(app.FileMenu,'Label','Close',...
                    'Separator','on','Callback',@app.closeApp);
            app.AnalysisMenu = uimenu(app.Figure,'Label','Analysis');
                uimenu(app.AnalysisMenu,'Label','Time Course Fitting',...
                    'Callback',@app.fitData);
                
            set([app.Figure,app.Axis,app.OpenBtn,app.FitBtn,...
                app.FilterBtn, app.FilterValTxt,app.AmpHistBtn,...
                app.BaseSubBtn,filterLabel,app.InvertChkBox],...
                'Units','normalized');
            
            app.NumFitPeriods = 0;
        end
        
        function closeApp(app,hObject,eventdata)
           selection = questdlg('End Model Fitting?',...
              'Exit Dialog',...
              'Yes','No','Yes'); 
           switch selection, 
              case 'Yes',
                 delete(app.Figure);
              case 'No'
              return 
           end
        end
        
        function openData(app,hObject,eventdata)
            if isempty(app.Filename)
                [file, path] = uigetfile({'*.abf','Axon binary file'});
            else
                path = fileparts(app.Filename);
                [file, path] = uigetfile({'*.abf','Axon binary file'},[],[path '\']);
            end
            if file == 0
                return
            end
            app.Filename = [path file];
            set(app.Figure,'Pointer','watch');
            drawnow;
            [app.Data,app.Header,app.Time] = abfread(app.Filename);
            plotdata(app);
            set(app.Figure,'Pointer','arrow');
            set(app.InvertChkBox,'value',get(app.InvertChkBox,'Min'));
        end
        
        function fitData(app,hObject,eventData)
            % Set the filter corner frequency (in sampling intervals)
            % for the step response
            user_entry = get(app.FilterValTxt,'string');
            fc = str2double(user_entry);
            if isnan(fc)
              msg=sprintf('"%s" is not a valid filter cutoff frequency\nYou must enter a numeric value',user_entry);
              errordlg(msg,'Bad Input','modal');
              uicontrol(app.FilterValTxt);
                return
            end
            adint = app.Header.fADCSequenceInterval*1e-3;
            fc = fc*adint; % fc = fc / (1/ts);
            app.Fc = fc;
            
            set(app.Figure,'Pointer','watch');
            drawnow;
            handles = guihandles(app.DisplayFigure);
            wtbar=waitbar(0,'Finding Transitions...','CreateCancelBtn',...
                'setappdata(gcbf,''canceling'',1)');
            setappdata(wtbar,'canceling',0);
            
            dat = app.Data;
            Tr = sqrt(log(2)/(2*pi)) / app.Fc;
            set(handles.riseTimeMsg,'string',num2str(Tr,3));
            
            [ start, fin, numOpenTrans, numShutTrans, openIdx, ...
                shutIdx, gm, I] = findtransitions (dat);
            baseline = max(gm.mu(:,1));
            dat = dat - baseline;
            openCurrent = min(gm.mu(:,1))-baseline;
            
            set(app.DisplayFigure,'Visible','on');
            
            amplitude = [];
            duration = [];
            begin = [];
            tbeg = [];
            tamp = [];
            tic;
            ini = 1;
            last = numel(start);
            for ii=ini:last
                % Check for Cancel button press
                if getappdata(wtbar,'canceling')
                    break
                end
                waitbar((ii-ini)/(last-ini),wtbar,sprintf('Fitting Open Period %d of %d',ii,last));
                
                dataToFit = dat(start(ii):fin(ii));
                numTrans = numel(openIdx{ii})+numel(shutIdx{ii});
                amps0 = zeros(numTrans,1);
                amps0(1:2:end) = openCurrent;
                trans0 = [openIdx{ii}-2, shutIdx{ii}-1]';
                trans0 = trans0(:);
%                 x0 = zeros(round(numTrans/2),3);
%                 t=[r(1,1); diff(r(:))];
%                 x0(:,1) = t(1:2:end);
%                 x0(:,2) = t(2:2:end);
%                 x0(:,3) = openCurrent;
                                
                lb_amps = 2*openCurrent*ones(size(amps0));
                lb_t = zeros(size(trans0));
                if ~isempty(begin)
                    lb_t(1) = begin(end) - start(ii);
                end
%                 ub_amps = 0.25*openCurrent*ones(size(amps0));
                ub_amps = -0.25*openCurrent*ones(size(amps0));
                ub_t = numel(dataToFit)*ones(size(trans0));
                % Fix the last transition to be closed (i.e. amplitude=0)
                idxFix = false(numTrans,1);
                idxFix(end) = true;
                try
                    [tamp, tbeg] = tcfixamp(dataToFit,amps0,trans0,...
                        app.Fc,find(idxFix),[lb_amps(~idxFix); lb_t],[ub_amps(~idxFix); ub_t],...
                        'plot',true,'figure',app.Figure);
                catch err
                    plot(app.Axis,dataToFit);
                    rethrow(err);
                    set(app.Figure,'Pointer','arrow');
                    delete(wtbar);
                    break;
                end
                % Test whether any transitions were less than 2x the filter
                % rise time, and if so fix those
                % but we don't know the duration of the last transition 
                % (because it is still ongoing) so do not fix it
                idxFix = [diff(tbeg)<2*Tr; false];
                
                OPENS = true(size(amps0));
                OPENS(2:2:end) = false;
                if any(idxFix)
                    idxGood = diff(tbeg)>=3*Tr;
                    if any(idxGood)
                        amps0(idxFix&OPENS) = mean(tamp(idxGood));
                        amps0(idxFix&~OPENS) = 0;
                    end
                    [tamp, tbeg] = tcfixamp(dataToFit, amps0, trans0,...
                        app.Fc, find(idxFix), ...
                        [lb_amps(~idxFix);lb_t],...
                        [ub_amps(~idxFix);ub_t],...
                        'plot',true,'figure',app.Figure,'showfit',false);
                end
%                 [rd,rs] = imposeres(diff(tbeg),tamp(1:end-1),0.053/0.025,0.032/0.025);
                amplitude = [amplitude; tamp];
                tbeg = tbeg + start(ii);
                tbeg(tbeg>fin(ii)) = fin(ii);
                begin = [begin; tbeg];
                if any(diff(begin)<0)
                    msgbox('Duration less than 0');
%                     disp(tbeg);
                    break;
                end
%                 waitbar((ii-ini)/(numel(start)-ini),wtbar);
                % Display results in the table
%                 try
%                     set(handles.fitResults,'Data',[tbeg-start(ii),diff([tbeg;tbeg(end)]),tamp]); %; ...
% %                         zeros(size(rd)), rd, rs]);
%                 catch
%                     msgbox('Displaying results did not work','Error');
%                     break
%                 end
%                 pause;
            end
            begin_t = adint*begin;
            duration = diff(begin_t);
            
            elapsedTime = toc;
            msg = sprintf('Fit %d transitions\nElapsed Time: %.1f s',ii,elapsedTime);
            msgbox (msg,'Finished');
            set(app.Figure,'Pointer','arrow');
            delete(wtbar);
            
            app.Durations = duration;
            app.Amplitudes = amplitude;
            app.TransTimes = begin_t;
            app.DwellStartIdx = start;
            app.DwellEndIdx = fin;
            app.OpenIdx = openIdx;
            app.ShutIdx = shutIdx;
            app.NumFitPeriods = ii;
        end
        
        function plotdata(app)
            axes(app.Axis);
            cla;
            title(app.Filename);
            plot(app.Time(1:1e2:end),app.Data(1:1e2:end));
            xlabel('Time (ms)');
            h=app.Header;
            ylabel([strtrim(char(h.sADCChannelName(1:10))), ' (', strtrim(char(h.sADCUnits(1:8))), ')']);
            [~, name] = fileparts(app.Filename);
            title(name,'Interpreter','none');
        end
        
        function filterData(app,hObj,eventdata)
            ON = get(hObj,'Value');
            if isempty(app.Data)
                set(hObj,'Value',~ON);
                errordlg('No Data','Error','modal');
                return
            end
            if ON
                user_entry = get(app.FilterValTxt,'string');
                Fc = str2double(user_entry);
                if isnan(Fc)
                  msg=sprintf('"%s" is not a valid entry\nYou must enter a numeric value',user_entry);
                  errordlg(msg,'Bad Input','modal');
                  uicontrol(app.FilterValTxt);
                  set(hObj,'Value',~ON);
                    return
                end
                Ts = app.Header.fADCSequenceInterval*1e-3;
                Fs = 1/Ts;
                Fc = Fc / Fs;
                app.Fc = Fc;
                sig = sqrt(log(2))/(2*pi*Fc);
                if sig < 0.62
                    numcoeffs = 1;
                else
                    numcoeffs = ceil(4*sig);
                end
                h_filt = gaussfir(Fc,numcoeffs,1);
                app.Data = filter(h_filt,1,app.Data);
                plotdata(app);
            else
                set(app.Figure,'Pointer','watch');
                drawnow;
                [app.Data] = abfread(app.Filename);
                plotdata(app);
                set(app.Figure,'Pointer','arrow');
            end
        end
        
        function filterVal(app,hObj,eventdata)
            user_entry = get(app.FilterValTxt,'string');
            numEntry = str2double(user_entry);
            if isnan(numEntry)
              msg=sprintf('"%s" is not a valid entry\nYou must enter a numeric value',user_entry);
              errordlg(msg,'Bad Input','modal');
              uicontrol(hObj);
              return
            end
        end
        
        function filterValBtnDwn(app,hObj,eventdata)
            selType = get(app.Figure,'SelectionType');
            if strcmpi(selType,'alt')
                msg = sprintf('Set filter corner frequency for applying a filter to the data\nor for use in time course fitting');
                msgbox(msg,'Help');
            end
        end
        
        function ampHist(app,hObj,eventData)
            if isempty(app.Data)
                errordlg('Open data first','No Data','modal');
                return
            end
            hFig = figure('Name','Amplitude Histogram','NumberTitle','off');
            hist(app.Data,min(app.Data):0.1:max(app.Data));
            h=app.Header;
            xlabel([strtrim(char(h.sADCChannelName(1:10))), ' (', strtrim(char(h.sADCUnits(1:8))), ')']);
            ylabel('Number');
        end
        
        function baselineCorrect(app,hObj,eventData)
            if isempty(app.Data)
                errordlg('Read data first', 'Error', 'modal');
                return
            end
            set(app.Figure,'Pointer','watch');
            drawnow;
            try
                idx = unique(randi(numel(app.Data),round(numel(app.Data)/100),1));
                dydt = diff([app.Data; app.Data(end)]);
                X = [app.Data(idx), dydt(idx)];
                gm = gmdistribution.fit(X,3);
                [~, state] = sort(gm.mu(:,1));
%                 app.Data = app.Data-gm.mu(state(3));
                I = cluster(gm,[app.Data, dydt]);
                closed = filter (ones(1,100)/100, 1, app.Data(I==state(3)));
                closed(1:99)=[];
                closed(end+1:end+99)=closed(end);
                baseline = interp1(app.Time(I==state(3)),closed,app.Time,...
                    'linear','extrap');
                app.Data = app.Data - baseline;
                axes(app.Axis);
                hold on;
                plot(app.Time,baseline,'c');
                hold off;
                pause(2);
                plotdata(app);
                set(app.Figure,'Pointer','arrow');
            catch err
                set(app.Figure,'Pointer','arrow');
                errordlg('Baseline Correct Did Not Work');
                rethrow(err);
            end
        end
        
        function saveIdl(app,hObj,eventData)
            if app.NumFitPeriods==0
                errordlg('No Idealized Data');
                return
            end
            
            if isempty(app.Filename)
                [file, path] = uiputfile({'*.dwt','QuB Idealized Data'});
            else
                [filename, path] = fileparts(app.Filename);
                filename = [filename, '.dwt'];
                [file, path] = uiputfile({'*.dwt','QuB Idealized Data'},[],[path filename]);
            end
            if file == 0
                return
            end
            
            adint = app.Header.fADCSequenceInterval*1e-3;
            dwtwrite(app.Durations,app.Amplitudes(1:end-1),adint,'filename',[path, file]);
        end
        
        function invertData(app,hObj,eventData)
            if ~isempty(app.Data)
                app.Data = -app.Data;
                plotdata(app);
            else
                set(hObj,'Value',get(hObj,'Min'));
            end
        end
        
        function exportToWksp(app,hObj,eventData)
            export2wsdlg({'Data','Time','Header'},...
                {'data','time','header'},...
                {app.Data,app.Time,app.Header});
        end
        
        function importFromWksp(app,hObj,eventData)
            vars = evalin('base','who');
            [datasel,ok] = listdlg('PromptString','Select data to import:',...
                'SelectionMode','single',...
                'ListString',vars);
            [timesel,ok] = listdlg('PromptString','Select time to import:',...
                'SelectionMode','single',...
                'ListString',vars);
            if isempty(datasel) || isempty(timesel)
                return
            end
            
            app.Data = evalin('base',vars{datasel});
            app.Time = evalin('base',vars{timesel});
            app.Header = [];
            app.Header.fADCSequenceInterval = (app.Time(2)-app.Time(1))*1e3;
            app.Header.sADCChannelName = 'Current';
            app.Header.sADCUnits = 'pA';
            plotdata(app);
        end
        
    end % methods
end