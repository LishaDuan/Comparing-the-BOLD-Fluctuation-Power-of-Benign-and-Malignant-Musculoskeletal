cd('E:\Tumor_frequency');   
load('E:\Tumor_frequency\TumorcomTC.mat');
SamplePeriod=2;
   for i=1:size(malignant,2);
       tumorTC(:,i)=malignant(:,i);
       
   %Plot the time course
% 	theTimeCourse =squeeze(AConfig.Dataset(theX, theY, theZ, :));
% 	axes(AConfig.hAxesTimeCourse);	cla;
	if rest_misc('GetMatlabVersion')>=7.3,	
		plot([1:length( tumorTC(:,i))] *SamplePeriod, ...
				tumorTC(:,i),'blue', 'DisplayName', 'Time Course');
	else	% Matlab 6.5 doesn't support  plot's property 'DisplayName'
		plot([1:length(tumorTC(:,i))] *SamplePeriod, ...
			tumorTC(:,i),'blue');
	end	
	xlim([1, length(tumorTC(:,i))] *SamplePeriod);
	theYLim =[min(tumorTC(:,i)), max(tumorTC(:,i))];
	if ~isfinite(theYLim(1)), theYLim(1)=0; end
	if ~isfinite(theYLim(2)), theYLim(2)=0; end
	if theYLim(2)>theYLim(1), ylim(theYLim); end	
	set(gca, 'Title',text('String','Time Course', 'Color', 'magenta'));
	xlabel('Time( seconds)');
	ylabel('Intensity');     
    
    
   %Plot the amplitude in Freq domain
	thePaddedLen(:,i) =rest_nextpow2_one35(length(tumorTC(:,i)));
				
		%Before FFT, remove linear trend first
		theTimeCourseNoTrend(:,i) =detrend(double(tumorTC(:,i)))  ;%...
					%+ repmat(mean(double(theTimeCourse)),[length(theTimeCourse), 1]);
                    
		%Draw the detrended data in the timeCourse Axes
% 		axes(AConfig.hAxesTimeCourse);	
%         hold on;
		thePlotTimeCourse(:,i) =theTimeCourseNoTrend(:,i) + repmat(mean(double(tumorTC(:,i))), [length(tumorTC(:,i)), 1]);
		if rest_misc('GetMatlabVersion')>=7.3,
			plot([1:length(tumorTC(:,i))] *SamplePeriod, ...
					thePlotTimeCourse(:,i), 'g:', 'DisplayName', 'Detrended Time Course');		
		else
			plot([1:length(tumorTC(:,i))] *SamplePeriod, ...
					thePlotTimeCourse(:,i), 'g:');		
		end
		theYLim =[min(theYLim(1),min(thePlotTimeCourse(:,i))), max(theYLim(2),max(thePlotTimeCourse(:,i)))];
		if ~isfinite(theYLim(1)), theYLim(1)=0; end
		if ~isfinite(theYLim(2)), theYLim(2)=0; end
		if theYLim(2)>theYLim(1), ylim(theYLim); end
		set(gca, 'Title',text('String','Time Course(Green dot line is after removing linear trend)', 'Color', 'magenta'));
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%Calculate the FFT
		thePowerTitle ='Power Spectrum after removing linear trend';
		
		theFreqSeries =fft(theTimeCourseNoTrend(:,i), thePaddedLen(:,i)); % multiply 2 just because I only plot half of the spectrum, so I make all enery to the plotted half part
					
	
	theSampleFreq		=1/SamplePeriod ;
	theFreqPrecision 	=theSampleFreq/thePaddedLen(:,i);
	theFreqLim =[theFreqPrecision: theFreqPrecision :theSampleFreq/2];
	theXLim =[2,(thePaddedLen(:,i)/2 +1)];	%don't Contain DC, because AFNI don't contain it in PowerSpectrum
	
	%Calcute the Power Spectrum
	theFreqSeries =abs(theFreqSeries([theXLim(1):theXLim(2)])); % Get the half's amplitude	
	theFreqSeries(1:end) =theFreqSeries(1:end).^2 /length(tumorTC(:,i));%Don't containt the DC component because abs didn't make DC 2-times to its original amplitude , dawnsong 20070629
	%theFreqSeries(1) =theFreqSeries(1) /length(theTimeCourse);	% now process the DC component
	
	%Since we dropped half the FFT, we multiply mx by 2 to keep the same energy.  
	% The DC component and Nyquist component, if it exists, are unique and should not 
	% be mulitplied by 2. 
	theFreqSeries(1:end-1) =theFreqSeries(1:end-1) *2;
    theFinalFreqSeries_magligant(:,i)=theFreqSeries;
  
	
% 	axes(AConfig.hAxesAmplitude);	
    cla;
	if rest_misc('GetMatlabVersion')>=7.3,
		plot(theFreqLim, theFreqSeries, ...
			'Color', 'blue', 'DisplayName', 'Power Spectrum');	
	else
		plot(theFreqLim, theFreqSeries, 'Color', 'blue');	
	end
	
	xlim([theFreqLim(1) , theFreqLim(end)]);
	theYLim =[min(theFreqSeries(1:end)), max(theFreqSeries(1:end))];
	if ~isfinite(theYLim(1)), theYLim(1)=0; end
	if ~isfinite(theYLim(2)), theYLim(2)=0; end
	if theYLim(2)>theYLim(1), ylim(theYLim); end	
	set(gca, 'Title',text('String',thePowerTitle, 'Color', 'magenta'));	
	xlabel(sprintf('Frequency( Hz), Sample Period( TR)=%g seconds', SamplePeriod));
	ylabel('Amplitude');
        
   end  
%    cd('E:\Tumor_frequency');   
%    save  theFinalFreqSeries_benighn  theFinalFreqSeries_benighn  
%    save  theFinalFreqSeries_magligant  theFinalFreqSeries_magligant 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AmplitudeSum_B=sum(theFinalFreqSeries_benigh);
    for ss=1:size(theFinalFreqSeries_benigh,2)
        for aa=1:size(theFinalFreqSeries_benigh,1)
            fabenign(aa,ss)=theFinalFreqSeries_benigh(aa,ss)./AmplitudeSum_B(1,ss);
        end
    end

     AmplitudeSum_M=sum(theFinalFreqSeries_magligant);
    for ss=1:size(theFinalFreqSeries_magligant,2)
        for aa=1:size(theFinalFreqSeries_magligant,1)
            famalignant(aa,ss)=theFinalFreqSeries_magligant(aa,ss)./AmplitudeSum_M(1,ss);
        end
    end
    save famplitude,  fabenign, famalignant 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% two-sample t-test
g_t=zeros(96,1);
g_p=zeros(96,1);
g_h=zeros(96,1);
for ff=1:96
     temp_b=fabenign(ff,:);
     temp_m=famalignant(ff,:);
    [h,p,ci,stat]=ttest2(temp_b,temp_m,0.05,'both');
        g_t(ff,1)=stat.tstat;
        g_p(ff,1)=p;
        g_h(ff,1)=h;
              
end
% save gc_p gc_p 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% two-sample t-test
gc_t=zeros(96,1);
gc_p=zeros(96,1);
gc_h=zeros(96,1);
for ff=1:96
     temp_b=theFinalFreqSeries_benigh(ff,:);
     temp_m=theFinalFreqSeries_magligant(ff,:);
    [h,p,ci,stat]=ttest2(temp_b,temp_m,0.05,'both');
        gc_t(ff,1)=stat.tstat;
        gc_p(ff,1)=p;
        gc_h(ff,1)=h;
              
end
% save gc_p gc_p 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonparamatic test: The Wilcoxon rank sum test (Mann-Whitney U test)
% gc_t=zeros(96,1);
n_p=zeros(96,1);
n_h=zeros(96,1);
for ff=1:96
    temp_b=fabenign(ff,:);
    temp_m=famalignant(ff,:);

%     [p h]=ranksum(temp_b, temp_m, 'alpha',0.05, 'method','exact');
      [p h]=ranksum(temp_b, temp_m, 'alpha',0.05);
      n_p(ff,1)=p;
      n_h(ff,1)=h;
end
save gc_p_ranksum gc_p 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot mean
% cd('E:\Tumor_frequency');
% load('theFinalFreqSeries_benigh.mat');
% load('theFinalFreqSeries_magligant.mat');
% meanfreBenign=mean(theFinalFreqSeries_benigh,2);
% meanfreMagligant=mean(theFinalFreqSeries_magligant,2);
% cla;
% theFreqLim =[theFreqPrecision: theFreqPrecision :theSampleFreq/2];
% plot(theFreqLim, meanfreBenign, ...
%     'Color', 'black', 'DisplayName', 'Power Spectrum');
% hold on
% plot(theFreqLim, meanfreMagligant, ...
%     'Color', 'blue', 'DisplayName', 'Power Spectrum');
% plot([0.1508 0.1508],[3,3],'r*');
% hold on
% plot([0.1534 0.1534],[3,3],'r*');
% hold on
% plot([0.247 0.247],[3,3],'r*');
% hold on
% 
% xlim([theFreqLim(1) , theFreqLim(end)]);
% theYLim =[min(meanfreBenign(1:end)), max(meanfreBenign(1:end))];
% if ~isfinite(theYLim(1)), theYLim(1)=0; end
% if ~isfinite(theYLim(2)), theYLim(2)=0; end
% if theYLim(2)>theYLim(1), ylim(theYLim); end
% set(gca, 'Title',text('String',thePowerTitle, 'Color', 'magenta'));
% % set(gca,'Marker','r*') 
% 
% xlabel(sprintf('Frequency( Hz), Sample Period( TR)=%g seconds', SamplePeriod));
% ylabel('Amplitude');





        
  
    
    
    

        
        
        
        
        
        
        