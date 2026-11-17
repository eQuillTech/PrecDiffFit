classdef Expt<handle
   properties
	name=                       '';
	filename=                   '';             
	H_refL=                     [];
    intensL=                    [];
    sigIL=                      [];
    FL=                         [];%Magnitudes of structure factors 
    plotL=                      [];%List of data peaks to plot
    doPlot=                     false;
    region;
    diff;
    selectData=                 true;
    lu_plotL=                   [];
    convEff=                    10;%Conversion efficiency of camera
    balanceIntens=              false;%use average of Friedel pairs
    includePairDiffError=       true;
    beamPairL=                  [];%Flag 1 if in pair
    fig;
    figNum;
    img;
   end
   
	methods
        function h=Expt(name)
            if nargin > 0
                h.name = name;
            end
        end
    
        %Find averages of Friedel pairs
        function expt=balance(expt)        
            nBeamExpt=size(expt.H_refL,1);
            expt.beamPairL=zeros(nBeamExpt,1);
            for ii=1:nBeamExpt
                H_refV1=expt.H_refL(ii,:);
                g1V=H_refV1*expt.region.B_refM;
                if(norm(g1V)<1e-7)
                    expt.beamPairL(ii,1)=1;
                    continue;
                end
                for jj=ii+1:nBeamExpt
                    H_refV2=expt.H_refL(jj,:);
                    g2V=H_refV2*expt.region.B_refM;
                    
                    gAvgV=g1V+g2V;
                    if norm(gAvgV)<1e-7 %%Remove duplicate peak
                        intens1=expt.intensL(ii,1);
                        intens2=expt.intensL(jj,1);
                        intensAvg=(intens1+intens2)/2;
                        if expt.balanceIntens
                            expt.intensL(ii,1)=intensAvg;
                            expt.intensL(jj,1)=intensAvg;
                            expt.sigIL(ii,1)=sqrt(intensAvg);
                            expt.sigIL(jj,1)=sqrt(intensAvg);
                       end
                        expt.beamPairL(ii,1)=1;
                        expt.beamPairL(jj,1)=1;
                       
                        %Include difference in error
                        if expt.includePairDiffError
                            sigI_Diff=abs(intens2-intens1)/2;
                            expt.sigIL(ii,1)=sqrt(sigI_Diff+expt.sigIL(ii,1)^2);
                            expt.sigIL(jj,1)=sqrt(sigI_Diff+expt.sigIL(jj,1)^2);
                        end
                    end
                end
            end
        end   
            
            
         %Read a text file containg a text hearder line followed by
        %four columns: h k l Intens
        function expt=readBeamL(expt,region)
            if expt.selectData
                [expt.filename,path] = uigetfile('*.txt');
                if isequal(expt.filename,0)
                   disp('File selection canceled');
                else
                   disp(['Opening ',fullfile(path,expt.filename)]);
                end
            end

            expt.H_refL=[];
            expt.intensL=[];
            expt.sigIL=[];

            % Read tab-delimited text (Octave-compatible), skipping the header row
            txt=fileread(expt.filename);
            lines=strsplit(strrep(txt,char(13),''),char(10));
            lines=lines(2:end);
            lines=lines(~cellfun(@(s) isempty(strtrim(s)),lines));
            nL=numel(lines);
            numericData=NaN(nL,4);
            for row=1:nL
                parts=strsplit(lines{row},char(9));
                for col=1:min(4,numel(parts))
                    numericData(row,col)=str2double(parts{col});
                end
            end
            expt.H_refL=numericData(:,1:3);
            expt.intensL=numericData(:,4);

          % HrefDataL=HrefDataL';
             %%Check peaks
            gL=expt.H_refL*region.B_refM;
            nBeamExpt=size(gL,1);
            beamFlag=ones(nBeamExpt,1);
            for ii=1:nBeamExpt
                gi=gL(ii,:);
                for jj=ii+1:nBeamExpt
                    dgV=gi-gL(jj,:);
                    if norm(dgV)<1e-7 %%Remove duplicate peaks
                        beamFlag(jj,1)=0;
                    end
                end
                if expt.intensL(ii)<0%remove negatives
                  beamFlag(ii)=0;
                end
            end

            H_refPL=zeros(0,3);
            intensPL=zeros(0,1);
            for ii=1:nBeamExpt
                if beamFlag(ii)==1
                    H_refPL=[H_refPL; expt.H_refL(ii,:)];
                    intensPL=[intensPL;expt.intensL(ii)];
                end
            end

            expt.H_refL=H_refPL;
            expt.intensL=intensPL;
            gL=expt.H_refL*region.B_refM;
            nBeamExpt=size(gL,1);

            %%Sort data peaks by length of g
            gL=expt.H_refL*region.B_refM;
            len_gL=vecnorm(gL')';
            [len_gL,sortL]=sort(len_gL);

            H_refPL=expt.H_refL(sortL(:),:);
            expt.H_refL=H_refPL;
            
            intensL=expt.intensL(sortL(:),:);           
            expt.intensL=intensL/expt.convEff;
            expt.sigIL=sqrt(expt.intensL);
            
            expt.balance();
           
            expt.plotL=(1:nBeamExpt)';%List of expt peaks to plot
        end
        
        %Create plot list
       function expt=createPlotList(expt,luL,plotFitBeamsOnly,gPlotMax)
            nBeam=size(expt.H_refL,1);
            expt.lu_plotL=[];
            nPlotBeam=0;
           for iBeam=1:nBeam
              H2_refV=expt.H_refL(iBeam,:)';                 
              g2V=expt.region.B_refM'*H2_refV;
              doPlotBeam=true;
             if norm(g2V)>gPlotMax
                 doPlotBeam=false;
             end
             if (plotFitBeamsOnly)&&(~ismember(iBeam,luL))
                  doPlotBeam=false;
               end
            if doPlotBeam
                nPlotBeam=nPlotBeam+1;
               expt.lu_plotL(nPlotBeam,1)=iBeam;
            end
          end
       end
       
         function expt=plot(expt,const)
            if isempty(expt.fig)
                expt.fig=figure('Name',expt.name,'NumberTitle','off');
                expt.figNum=get(gcf,'Number');

            end
  
            if expt.diff.doRealPlot
                expt.img=expt.diff.realPlot(const,expt.region,expt.intensL,expt.H_refL,expt.lu_plotL);
            else
               expt.img=expt.diff.plot(const,expt.region,expt.intensL,expt.H_refL,expt.lu_plotL);
            end   
                showFig=false;
            if ishghandle(gcf)
                frontFigNum=get(gcf,'Number');
                if frontFigNum~=expt.figNum
                    showFig=true;
                end
            end
            if showFig
                figure(expt.figNum);
            end
            imshow(expt.img,'InitialMagnification',100);     
        end


    end

end
