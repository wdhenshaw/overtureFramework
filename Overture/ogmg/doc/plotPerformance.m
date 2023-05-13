
%
% Plot performance comparisons
%
% Examples:
%   plotPerformance -case=solve -verticalBars=0 -savePlots=1
% 
%   plotPerformance -case=memory -verticalBars=1 -relative=0
%   plotPerformance -case=memory -verticalBars=0 -relative=0 -savePlots=1
%
% Initial version: Jul 9, 2022
%
% cd ~/Dropbox/CG6backup/overtureFramework/Overture/ogmg/doc
% cd ~/Dropbox/cg66backup/overtureFramework/Overture/ogmg/doc
%

function plotPerformance(varargin)

caseOption = 'solve';
savePlots=0;
verticalBars=0; % vertical or horizontal bars
yMax=0; yMin=1;
relative=1;     % show relative reults, normalized by Ogmg

% --- read command line args ---
for i = 1 : nargin
  line = varargin{i};
  caseOption       = getString( line,'-caseOption',caseOption );
  caseOption       = getString( line,'-case',caseOption );

  relative         = getInt( line,'-relative',relative );
  savePlots        = getInt( line,'-savePlots',savePlots );
  verticalBars     = getInt( line,'-verticalBars',verticalBars );
  yMax             = getReal( line,'-yMax',yMax );
  yMin             = getReal( line,'-yMin',yMin );
end 



% Set defaults for plotting 
fontSize=24; lineWidth=2; markerSize=8; 
set(0,'DefaultLineMarkerSize',markerSize);
set(0,'DefaultLineLineWidth',lineWidth);
set(0,'DefaultAxesFontSize',fontSize);
set(0,'DefaultLegendFontSize',fontSize);





% --- LOAD PERFORMANCE DATA -----

perfData;

% data


% relative CPU times
% data(cpuTotal,  2  ,1,square1024)=    0.58  ; % data(:,order,solver,gridType)= value

gridName{1}="square"; 
gridName{2}="cic";    
gridName{3}="shapes"; 

% quant=cpuTotal; 
if( strcmp(caseOption,'solve')==1 )
  quant=cpuSolve; quantName{1}="solve"; axisLabel='CPU time'; plotName='relativeSolveTimes'; 
  myTitle=sprintf('Relative solve times');
  gMax(1)=150; % 180; 
  gMax(2)=140; 
  gMax(3)=60; 
elseif( strcmp(caseOption,'memory')==1 )
  quant=storage; quantName{1}="memory";  
  if relative==1 
    myTitle=sprintf('Relative memory usage'); axisLabel='Memory'; plotName='relativeMemoryUsage';
  else
    myTitle=sprintf('Memory usage'); axisLabel='Memory reals/point'; plotName='memoryUsage';
  end 
  gMax(1)=200; % 25; 
  gMax(2)=200; % 20; 
  gMax(3)=200; % 15;   
else
  fprintf('ERROR: unknown case=%s\n',caseOption); pause; pause;
end

for gridType=1:3 % grids
  figure(gridType);

  for group=1:2
    ord = 2*group;
    if relative==1 
      cpu(group,1:4) = data(quant,ord,1:4,gridType)./data(quant,ord,1,gridType);
    else
      cpu(group,1:4) = data(quant,ord,1:4,gridType);
    end
  end

  % if 1==1 || verticalBars
  %   group=1; ord=2; 
  %   cpu(group,1:4) = data(quant,ord,1:4,gridType)./data(quant,ord,1,gridType);
  %   group=2; ord=4;
  %   cpu(group,1:4) = data(quant,ord,1:4,gridType)./data(quant,ord,1,gridType);; 
  % else
  %   group=1; ord=2; 
  %   cpu(group,4:-1:1) = data(quant,ord,1:4,gridType)./data(quant,ord,1,gridType);
  %   group=2; ord=4;
  %   cpu(group,4:-1:1) = data(quant,ord,1:4,gridType)./data(quant,ord,1,gridType);; 
  % end

  
  if verticalBars 
    hb = bar(cpu(1:2,1:4));
  else
    hb = barh(cpu(1:2,1:4));
  end
  % hb(1).FaceColor = [.2 .6 .5];
  hb(1).FaceColor = [.2 .6 .5]*1.4; % green 
  hb(2).FaceColor = [1 .5 .5];
  hb(3).FaceColor = [.5 .6 1];
  % hb(4).FaceColor = [.3 .2 .7]*1.4; % purple
  hb(4).FaceColor = [.5 .4 1]; % purple
  title(sprintf('%s, %s',myTitle,gridName{gridType}));
  str = {'Order 2'; 'Order 4'};
  if verticalBars 
    set(gca, 'XTickLabel',str, 'XTick',1:numel(str));
    ylabel(axisLabel);
  else
    set(gca, 'YTickLabel',str, 'YTick',1:numel(str)); ytickangle(90);
    xlabel(axisLabel);
  end
  legend('Ogmg','AMG','BiCGSt','GMRES','Location','best');
  grid on; 
  if verticalBars==0 
    yMax=gMax(gridType); 
    if yMax>0 
      xlim([0,yMax]);
    end
  else
    yMax=gMax(gridType); 
    if yMax>0 
      ylim([0,yMax]);
    end
  end

    % -------------------- LABEL TOP OF BARS -------------
    relCpu = cpu(1:2,1:4);
    xShift0=.15; 
    labelSize=20;
  
    labelBars( hb,relCpu,yMax,xShift0,labelSize,verticalBars );

  if savePlots 
    fprintf('Pause before saving plot to adjust legend...\n');
    pause
    if verticalBars 
      fullPlotName=sprintf('fig/%s%sV2',plotName,gridName{gridType});
    else
      fullPlotName=sprintf('fig/%sHorizontalBars%sV2',plotName,gridName{gridType});
    end
    savePlotFile(fullPlotName,'pdf');   
  end
end % grids 

return
end


% ----------------------------------------------------
% Function to put numbers on top of bar chart bars
% ----------------------------------------------------
function labelBars( hb,relCpu,yMax,xShift0,labelSize,verticalBars )

  nd=2; 

  M=relCpu;
  numbersToAdd=relCpu;
  barWidth = hb.BarWidth;
  numCol = size(M,2);
  cnt = 0;
  for ii = numbersToAdd'
      cnt = cnt + 1;
      xPos = linspace(cnt - barWidth/2, cnt + barWidth / 2, numCol+1);
      idx = 1;
      group=0;
      for jj = xPos(1:end-1)
          group=group+1; 
          val = numbersToAdd(cnt,idx);
          y = min(yMax-15,M(cnt,idx));  % restrict position  if the label goes off plot
          % xShift=(.08 - (idx-1)*.02)*barWidth; 
          xShift=(xShift0 - (idx-1)*.02)*barWidth; 
          % xShift=(.275 - (idx-1)*.085)*barWidth; 
          if( nd ==2 ) yShift=2; else yShift=2.7; end;  
          % fprintf('idx=%g, jj=%g xShift=%g num=%s\n',idx,jj,xShift,num2str(val,'%0.1f'));
          if verticalBars
            if val==1 || group==1
              ht=text(jj+xShift, y + yShift, sprintf('%0.0f Ogmg (best)',val),'fontSize',labelSize); set(ht,'Rotation',90);
            else
              ht=text(jj+xShift, y + yShift, num2str(val,'%0.0f'),'fontSize',labelSize); set(ht,'Rotation',90);
            end
          else
            % horizontal bars
            if val==1 || group==1
              ht=text(y + yShift, jj+xShift, sprintf('%0.0f Ogmg (best)',val), 'fontSize',labelSize); set(ht,'Rotation',0);
            else
              ht=text(y + yShift, jj+xShift,  num2str(val,'%0.0f'),'fontSize',labelSize); set(ht,'Rotation',0);
            end            
          end
          idx = idx +1;
      end     
  end
  return
end

% if 1==0 
%   % **** OLD ***
% figure(1);


% % y = [2 2 3; 2 5 6; 2 8 9; 2 11 12];
% cpu = [ .58, 12.8; .82, 13.8; 1.36, 25.6; 1.9, 112]; % cou times for Ogmg, AMG; ...


% str = {'sq O2'; 'sq O4'; 'cic O2'; 'cic O4'};
% hb = bar(cpu);
% legend('Ogmg','AMG','Location','north');
% title('CPU: Ogmg and AMG')
% set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
% grid on; ylabel('seconds');
% yMax=120; 
% ylim([0,yMax]);

%  % -------------------- LABEL TOP OF BARS -------------
%   relCpu = cpu;
%   nd=2; 

%   M=relCpu;
%   numbersToAdd=relCpu;
%   barWidth = hb.BarWidth;
%   numCol = size(M,2);
%   cnt = 0;
%   for ii = numbersToAdd'
%       cnt = cnt + 1;
%       xPos = linspace(cnt - barWidth/2, cnt + barWidth / 2, numCol+1);
%       idx = 1;
%       for jj = xPos(1:end-1)
%           val = numbersToAdd(cnt,idx);
%           y = min(yMax-10,M(cnt,idx));
%           % xShift=(.08 - (idx-1)*.02)*barWidth; 
%           xShift=(.275 - (idx-1)*.085)*barWidth; 
%           if( nd ==2 ) yShift=2; else yShift=2.7; end;  
%           % fprintf('idx=%g, jj=%g xShift=%g num=%s\n',idx,jj,xShift,num2str(val,'%0.1f'));
    
%           labelSize=16;
%           ht=text(jj+xShift, y + yShift, num2str(val,'%0.1f'),'fontSize',labelSize); set(ht,'Rotation',90);
%           idx = idx +1;
%       end     
%   end

% % ----------------
% figure(2)
% cpu = [ .58/12.8; .82/13.8; 1.36/25.6; 1.9/112]; % cou times for Ogmg, AMG; ...
% cpu = 1./cpu; 

% str = {'sq O2'; 'sq O4'; 'cic O2'; 'cic O4'};
% hb = bar(cpu);
% legend('AMG/Ogmg','Location','north');
% title('Ogmg CPU speedup factor vs AMG')
% set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
% grid on; ylabel('AMG/Ogmg');


% end


% if 1==0 
%     relCpu = [totalCpuArray(3,:)./totalCpuArray(1,:); totalCpuArray(4,:)./totalCpuArray(2,:) ];
    
%     % bar(diag(relCpu),'stacked');
%     hb=bar(relCpu);
%     grid on;
%     title(sprintf(' Rel-Cost : Curvilinear/Cartesian, %s',myTitle));
%     if( numCases ==4 )
%       legend('FD','FDA','UPC','UW','Location','north');
%       str = {'Order 2'; 'Order 4'};
%     else
%       legend('FD','FDA','UPC','UW','Location','northwest');
%       str = {'Cart/Curv'};
%     end 
%     set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
%     set(gca,'fontSize',fontSize);
    
%     % -------------------- LABEL TOP OF BARS -------------
%     M=relCpu;
%     numbersToAdd=relCpu;
%     barWidth = hb.BarWidth;
%     numCol = size(M,2);
%     cnt = 0;
%     for ii = numbersToAdd'
%         cnt = cnt + 1;
%         xPos = linspace(cnt - barWidth/2, cnt + barWidth / 2, numCol+1);
%         idx = 1;
%         for jj = xPos(1:end-1)
%             val = numbersToAdd(cnt,idx);
%             y = M(cnt,idx);
%             xShift=(.08 - (idx-1)*.02)*barWidth; 
%             if( nd ==2 ) yShift=2; else yShift=2.7; end;  
%             % fprintf('idx=%g, jj=%g xShift=%g num=%s\n',idx,jj,xShift,num2str(val,'%0.1f'));
      
%             text(jj+xShift, y + yShift, num2str(val,'%0.1f'),'fontSize',fontSize);
%             idx = idx +1;
%         end     
%     end
    
    
%     plotFileName=sprintf('costCurvilinearToCartesian%s.eps',caseLabel);
%     print('-depsc',plotFileName);
%     fprintf('Saved plot %s\n',plotFileName);

%     commandwindow; pause;
% end