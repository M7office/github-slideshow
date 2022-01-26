%%%%%%%%%%%%%%%%%%%%%%%%
%figure left is the whole analysis; right part focuses on specific voltage
% need to have same win start point and interval with
% cavitationAnalysisCathScript.m
% doubelcheck VolData
% 1 pulse or several pusles in each exp--> check the datapreprocess
% function
% check the zero padding
% awgFreq Mhz other Khz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
load('combined5030khz.mat');%li's
[~,index] 	= sortrows([voltageData.voltage].'); voltageData = voltageData(index); clear index

caseNum = 0;%li's
winNum  = 0;%li's
vol     = 120;%li's 


awgFreq    = voltageData(1).frequency/1000;   
delayTime  = voltageData(1).delayTime ;
sampleFreq = voltageData(1).sampleFreq;

	if (caseNum == 1) || (caseNum == 2) || (caseNum == 3) ||(caseNum == 4)
		cycleTime  = 10; %time for 50 cylces in us
	else if (caseNum == 5) || (caseNum == 6) || (caseNum == 7) ||(caseNum ==8)
		cycleTime  = 100;
	else if caseNum ==9
		cycleTime  = 1000;
	else if caseNum==0
		cycleTime = (1/awgFreq)*50;
		end
		end
		end
	end
	
	if awgFreq 			== 0.68 
		startPoint      = (-delayTime*sampleFreq)+200;%li's
		predictFactor   = 5.972;
		ICdoseStartFreq = 9; %li's 
		ICdoseEndFreq   = 10; %li's 
		cycleTime       = (1/awgFreq)*7;%li's % 10us 50cyc for 5030Mhz; 10us 20cyc for 1915khz;
	else if awgFreq 	== 1.915 
		startPoint      = (-delayTime*sampleFreq)+500;%li's
		predictFactor   = 10.4;
		ICdoseStartFreq = 12; %li's 
		ICdoseEndFreq   = 13; %li's 
		cycleTime       = (1/awgFreq)*20;%li's % 10us 50cyc for 5030Mhz; 10us 20cyc for 1915khz;
	else if awgFreq 	== 5.03 
		startPoint		= (-delayTime*sampleFreq)+1000;%li's
		predictFactor	= 9.43;
		ICdoseStartFreq = 18; %li's 
		ICdoseEndFreq 	= 19; %li's 
	else if awgFreq 	== 0.377
		predictFactor 	= 9.404;%li's
	else if awgFreq 	==4.35
		predictFactor 	= 31.956;
	else if awgFreq 	== 4.96
		predictFactor 	= 17.73;
	else if awgFreq 	== 5.565
		predictFactor 	= 27.038;
	else if awgFreq 	== 0.705
		predictFactor 	= 6;
	else if awgFreq 	== 1.985
		predictFactor 	= 9.747;
		end
		end
		end
		end
		end
		end
		end
		end
	end

winInterval 	= cycleTime*sampleFreq;
endPoint 	 	= startPoint+  winInterval;     
endPointNoise   = length(voltageData(1).waveforms)-1;
startPointNoise = endPointNoise- winInterval;




	if isfield(voltageData,'waveforms')
		new = 1;
	else
		new = 0;
	end


warning('off','MATLAB:colon:nonIntegerIndex')

	if new
		clear voltageDataNew
		k=1;
		for i=1:length(voltageData)
			pulseNum = size(voltageData(i).waveforms,2);
				for j=1:pulseNum
					voltageDataNew(k) 			= voltageData(i);
					voltageDataNew(k).waveforms = voltageData(i).waveforms(:,j);
					k 							= k+1;
				end
		end
		for kk = 1:length(voltageDataNew)
			voltageDataNew(kk).waveform 		= voltageDataNew(kk).waveforms;
		end
	    voltageData = voltageDataNew;
	end



hydroFreqList 			= 250:50:12000;
% hydroDataFP188_15P	= [0.284641545	0.196588775	0.118245179	0.108179969	0.099444617	0.104044973	0.095267226	0.09195307	0.090992967	0.090829352	0.092666855	0.087401612	0.092235592	0.089170937	0.096937691	0.095977438	0.093225268	0.097576626	0.096732963	0.096741256	0.098714669	0.099887731	0.101769669	0.103408519	0.103397797	0.104454837	0.10721211	0.10585945	0.108366204	0.111979801	0.113402664	0.116161079	0.116098579	0.116047839	0.118646963	0.119573953	0.124015459	0.121469345	0.126268562	0.128736096	0.128017474	0.127562731	0.129974845	0.1335907	0.131862816	0.129669032	0.136059741	0.135677296	0.135880976	0.139683311	0.138485953	0.141182917	0.140483047	0.141055139	0.141682945	0.145666288	0.117087467	0.116937022	0.35153673	0.608728826	0.461252513	0.15381532	0.165849805	0.149859254	0.16075054	0.159896145	0.168267728	0.17390341	0.178466146	0.180242324	0.178966776	0.187278669	0.192558521	0.188822278	0.191326138	0.203083077	0.206931987	0.204685712	0.206532568	0.217990034	0.221847	0.223098742	0.222145147	0.22006188	0.22279755	0.223648937	0.227542651	0.231135873	0.231247907	0.23143829	0.23624398	0.241902045	0.241383232	0.246115666	0.249888413	0.251677772	0.254904995	0.25850309	0.264295445	0.263595152	0.26960515	0.277157158	0.276466098	0.270712033	0.274380181	0.285778453	0.287606061	0.288313323	0.291308729	0.295630793	0.296238921	0.295334058	0.304721488	0.300665777	0.302016843	0.314364019	0.323582374	0.326792361	0.324030848	0.326930516	0.326500486	0.330230434	0.334026599	0.337535113	0.338228769	0.332708571	0.331828262	0.333361177	0.336433644	0.346611705	0.345416952	0.348694215	0.346968222	0.347428273	0.349117732	0.352673841	0.352772436	0.355771934	0.362272617	0.367563277	0.366890683	0.368587304	0.361478421	0.362146712	0.359883382	0.366229298	0.374570049	0.373829847	0.370975914	0.36544873	0.370455778	0.366123554	0.371719712	0.370453772	0.375442824	0.378083096	0.380265	0.383088298	0.381924337	0.385974981	0.390935087	0.392805654	0.387394372	0.386728281	0.38045529	0.37794832	0.383838434	0.381118176	0.377249874	0.382461808	0.386693268	0.381990642	0.386237244	0.3891189	0.384326331	0.390795739	0.392827063	0.386184832	0.387457704	0.384155995	0.382588614	0.388890057	0.385914464	0.381629856	0.385835078	0.382509797	0.37302943	0.374027952	0.381403754	0.389925823	0.386322728	0.392558061	0.387706248	0.391761773	0.392243529	0.395785883	0.405402865	0.410703683	0.389928238	0.403868187	0.398291288	0.39794009	0.401730719	0.396913611	0.402874304	0.393410645	0.402865473	0.399520595	0.395887918	0.399454499	0.409447333	0.399828579	0.39198527	0.381890882	0.39539692	0.390408637	0.392545019	0.388521333	0.390895844	0.395496157	0.400738389	0.395564522	0.410802448	0.413539246	0.413382296	0.42599629	0.406957832	0.416026163	0.402881917	0.405265151	0.425823292	0.393182923	0.43514173	0.409401055	0.439649715	0.42776003];
%0712%li's
hydroDataFP188_14P		= [0.149895875107240;0.137960987346350;0.0968073242745783;0.108579526875482;0.110382252872128;0.100337068887595;0.107490395530857;0.0976094966298150;0.112369914204269;0.106404322278313;0.103214112861148;0.103592841456619;0.105416479996660;0.105831433842281;0.107642178421141;0.112351717485408;0.104376888563304;0.114397103562942;0.108427714701800;0.110789572282186;0.111441068149510;0.110507877326176;0.116660846463431;0.115977355249387;0.117075599975629;0.119561182833191;0.119776099138940;0.121716825236071;0.123232251653171;0.124515100150018;0.128412509109005;0.127049256330662;0.129889576458841;0.132699240877677;0.131694881585321;0.132300311000298;0.136924165623989;0.135258444015289;0.140236594693967;0.142013550243616;0.142743650602491;0.144339564435007;0.144868772926729;0.147949164671521;0.148544800397664;0.147038874833119;0.150590987701657;0.151632485813129;0.152783851736893;0.155028218037921;0.155893603954794;0.159453923058453;0.159938560909247;0.162088990756867;0.163772391464836;0.165860052378505;0.165763909668897;0.169016669757295;0.224753256135197;0.227470342898016;0.104437654325513;0.162899750243387;0.180340998596655;0.184292774725279;0.182905362768123;0.185725845360037;0.190639856131219;0.194652549547879;0.197233877540432;0.199692856327341;0.199981270333952;0.204695382712644;0.205579944555796;0.205337993952769;0.212715975499559;0.216377075151644;0.220993072009176;0.221882846245293;0.223967324732790;0.229802347897685;0.234513764641494;0.237350630851369;0.236480433563232;0.235998650736546;0.239169519337874;0.242857712875355;0.243487362569530;0.241863162233374;0.246950856919087;0.249016125367951;0.253922756871411;0.258301255444974;0.260888302630023;0.264071175513943;0.265983529069774;0.270558182555932;0.273438394084421;0.274663797208037;0.278260167658507;0.280303169192428;0.284230272509812;0.292356998118647;0.293662634367372;0.290419431157750;0.293741146584959;0.296034589868415;0.300675641748467;0.309588494975317;0.312237831834898;0.313906971945244;0.320177142871000;0.316440500184465;0.316375109197906;0.316489841693891;0.319813743604068;0.324721835522416;0.327810260511560;0.334662318399651;0.340427394678324;0.340425173000709;0.341195643841663;0.343432731363104;0.347842675061463;0.353442375661371;0.354003950096181;0.351599348434000;0.349240499076085;0.350625426988628;0.353665376287630;0.358681006653426;0.362354631559916;0.364593205902462;0.364699791073779;0.368503407911156;0.368721178737980;0.370451327086694;0.369703417138982;0.368867588749092;0.370203712264220;0.379110661919363;0.372337476923093;0.375628673733973;0.374358564444605;0.371541668456415;0.371667563058565;0.374862694303360;0.377172976275748;0.374986983887105;0.376537516122023;0.375081161809636;0.377900509600687;0.372354633049307;0.376670256597466;0.380596173606559;0.378097504970076;0.374238241954926;0.381172324573296;0.380216565594480;0.382908269699613;0.380224391574554;0.378852008898875;0.380441208334206;0.372968562555790;0.378609003984536;0.378299012083013;0.379670043367071;0.377799061788781;0.375039134400030;0.375987933073374;0.364845693008432;0.369823881967264;0.364891069827061;0.360351899265243;0.367969249432114;0.361095265641634;0.368696881552344;0.365625430762693;0.363624173080055;0.364648618609107;0.364578419909084;0.358266468986950;0.361991389684785;0.353083191447825;0.347685029885005;0.350383887826597;0.335201180378467;0.321317941580889;0.346588701047185;0.350061162300587;0.345898782601523;0.341337943320756;0.345884432017794;0.341513652020745;0.339512290250810;0.337173020812053;0.341469905840701;0.342510473778564;0.340854786479628;0.337450430653776;0.330320683023582;0.329661473763793;0.333668238792664;0.328663878399557;0.319106973319115;0.321509794148613;0.312464508923658;0.317045195697518;0.311738700600173;0.310269208943179;0.311709952405998;0.310008733334540;0.299522554906643;0.303777862605082;0.304318165236708;0.301167370621885;0.302200746644193;0.306623942928187;0.305471266363199;0.303369117141776;0.308789131329335;0.319029682906421;0.314293025886402;0.320644773402126;0.333886654709630;0.336132045588293;0.344694406647861;0.348816083747079;0.347715565344421;0.352406098911894;0.357446576448931;0.371305681689081;0.383206142347925;0.381627823984274;0.382090254313049;0.384434235310375;0.394867606019812];
%0712%li's
closestFreq 			= max(round(awgFreq*1000/50)*50,250); %nearest 50kHz
factor 					= hydroDataFP188_14P(find(hydroFreqList == closestFreq,1));

	for i=1:length(voltageData)   
		fourierData							= voltageData(i).waveform(startPoint:endPoint);
		fourierDataNoise 					= voltageData(i).waveform(startPointNoise:endPointNoise);
		[voltageData(i).fourierData] 		= fourierData;
		[voltageData(i).fourierDataNoise] 	= fourierDataNoise;
		[voltageData(i).P1,voltageData(i).f]= computeFFT(hanning(length(fourierData)).*fourierData, voltageData(i).sampleFreq); %
		[voltageData(i).P1noise, ~] 		= computeFFT(hanning(length(fourierDataNoise)).*fourierDataNoise, voltageData(i).sampleFreq); %   
	end
	


clear IC_dose_mean_dBnoise_std mean_meanPressure std_meanPressure meanFFT IC_dose_mean_dB IC_dose_mean meanFFTnoise IC_dose_mean_dBnoise IC_dose_mean_dB_std
voltList = unique(nonzeros([voltageData.voltage]));

	for volIndex = 1:length(voltList)
	
		numTimes(volIndex)  				= size(find([voltageData.voltage] == voltList(volIndex)),2);
		startPointIdx       				= find([voltageData.voltage] == voltList(volIndex),1,'first');
		endPointindex   					= find([voltageData.voltage] == voltList(volIndex),1,'last');

		meanFFT(volIndex,:) 				= mean([voltageData(startPointIdx:endPointindex ).P1],2);
		meanFFTnoise(volIndex,:)			= mean([voltageData(startPointIdx:endPointindex ).P1noise],2); 
		[~,fStartIndex] 					= min(abs(voltageData(1).f-ICdoseStartFreq));
		[~,fEndIndex]   					= min(abs(voltageData(1).f-ICdoseEndFreq));
	
		IC_dose_mean_dB(volIndex)          	= mag2db(mean(meanFFT(volIndex,fStartIndex:fEndIndex))); %149 = 4.0224 MHz, 517 = 14.024 MHz
		IC_dose_mean_dB_std(volIndex)      	= std(mag2db(meanFFT(volIndex,fStartIndex:fEndIndex))); %149 = 4.0224 MHz, 517 = 14.024 MHz
		IC_dose_mean_dBnoise(volIndex)     	= mag2db(mean(meanFFTnoise(volIndex,fStartIndex:fEndIndex))); %149 = 4.0224 MHz, 517 = 14.024 MHz
		IC_dose_mean_dBnoise_std(volIndex) 	= std(mag2db(meanFFTnoise(volIndex,fStartIndex:fEndIndex))); %149 = 4.0224 MHz, 517 = 14.024 MHz    


		sampleData 							= voltageData(startPointIdx:endPointindex); 
		[press_yline,press_y800]         	= DataPreProcess(sampleData);
		[pressureY800(volIndex,:)] 			= (press_y800/factor)*predictFactor/1000;
		[PressureYline(volIndex,:)] 		= (press_yline/factor)/1000;
		mean_Pressure(volIndex)     		= mean(pressureY800(volIndex,:));
		std_Pressure (volIndex)     		= std (pressureY800(volIndex,:));
	end


fullfig;
subplot(4,2,1)
yyaxis left
errorbar(voltList,IC_dose_mean_dB,IC_dose_mean_dB_std)

%determine threshold
[~, maxDiffIdx] = max(diff(IC_dose_mean_dB));
hold on
scatter(voltList(maxDiffIdx),IC_dose_mean_dB(maxDiffIdx))
scatter(voltList(maxDiffIdx+1),IC_dose_mean_dB(maxDiffIdx+1));
scatter(voltList(find(voltList==vol)),IC_dose_mean_dB(find(voltList==vol)));

hcurPoint = scatter(voltList(end),IC_dose_mean_dB(end),'*');
% text(voltList(end),IC_dose_mean_dB(end),num2str(round(IC_dose_mean_dB(end)-IC_dose_mean_dBnoise(end))));
text(voltList(maxDiffIdx+1),IC_dose_mean_dB(maxDiffIdx+1),num2str(round(IC_dose_mean_dB(maxDiffIdx+1)-IC_dose_mean_dBnoise(maxDiffIdx+1))));
text(voltList(find(voltList==vol)),IC_dose_mean_dB(find(voltList==vol)),num2str(round(IC_dose_mean_dB(find(voltList==vol))-IC_dose_mean_dBnoise(find(voltList==vol)))));


errorbar(voltList,IC_dose_mean_dBnoise,IC_dose_mean_dBnoise_std)
title([num2str(voltageData(1).frequency) ' KHz'])
xlabel('Applied voltage [V]')
ylabel('IC Dose')


yyaxis right
errorbar(voltList,mean_Pressure,std_Pressure);
ylabel('Predicted Pressure [Mpa]')
grid on
grid minor
% for jj = 1:length(voltList)
% text(voltList(jj),mean_Pressure(jj),num2str(round(mean_Pressure(jj),1)),'FontSize',6);
% end
% text(voltList(end),mean_Pressure(end),num2str(round(mean_Pressure(end),1)));
text(voltList(maxDiffIdx+1),mean_Pressure(maxDiffIdx+1),num2str(round(mean_Pressure(maxDiffIdx+1),1)));
scatter(voltList(maxDiffIdx+1),mean_Pressure(maxDiffIdx+1),'*');
% scatter(voltList(end),mean_Pressure(end),'*');
text(voltList(find(voltList==vol)),mean_Pressure(find(voltList==vol)),num2str(round(mean_Pressure(find(voltList==vol)),1)));%li's
scatter(voltList(find(voltList==vol)),mean_Pressure(find(voltList==vol)),'*');%li's





subplot(4,2,3)
hold on
plot(voltageData(1).f,mag2db(meanFFT(maxDiffIdx,:)))
plot(voltageData(1).f,mag2db(meanFFT(maxDiffIdx+1,:)))
xlim([0 20])
ylim([-100 -20])
line([ICdoseStartFreq ICdoseStartFreq],[mean(IC_dose_mean_dBnoise)*.9 mean(IC_dose_mean_dBnoise)*1.1],'Color','g');
line([ICdoseEndFreq ICdoseEndFreq],[mean(IC_dose_mean_dBnoise)*.9 mean(IC_dose_mean_dBnoise)*1.1],'Color','g');
xlabel('Frequency [MHz]')
ylabel('dB')


subplot(4,2,5)
hold on

	for k=1:size(meanFFT,1)
		plot(voltageData(1).f,mag2db(meanFFT(k,:)))
	end
	
xlim([0 20])
ylim([-100 -20])
xlabel('Frequency [MHz]')
ylabel('dB')
title('FFTs of different volt');




subplot(4,2,7)
YList= unique([voltageData.y]);
minY = YList(1);
maxY = YList(end); 
	if length(YList)==1
		YStep = 1;    
	else    
		YStep = YList(2)-YList(1);
	end

x = (minY:YStep:maxY)/400;
	for i = 1:length(voltList)	
		scatter(x,PressureYline(i,:));
		hold on
	end
	
ylabel('Pressure [Kpa]');
xlabel('Axial Distance [mm]');
title('Yline Scan');
%ylim([0 4000])
xlim([-3.5 4])
hold off


startindex 	= find([voltageData.voltage] == vol,1,'first');
endindex 	= find([voltageData.voltage] == vol,1,'last');
VolData 	= voltageData(startindex:endindex);

	for i=0:winNum
		runAnalysis(VolData,i,startPoint,winInterval/(winNum+1),i);
		% pause(2)
	end
	
pulseByPulse(VolData,winNum,startPoint,winInterval);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(4,2,6)
hplot = plot(voltageData(end).waveform);%hold on;plot(voltageData(maxDiffIdx+1).waveform)
hold on
xline(startPoint);
xline(endPoint);
xline(startPointNoise);
xline(endPointNoise);
ylim([-0.8 0.8]);
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(4,2,8)
hplot2 = plot(voltageData(end).f,voltageData(end).P1);
xlim([0 20])
ylim([-100 -20])

% b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],'value',1, 'min',1, 'max',length(voltageData));
b = uicontrol('Style','slider','Position',[81,54,419,23],'value',1, 'min',1, 'max',length(voltageData));
b.SliderStep = [1/length(voltageData), 10/length(voltageData)];
b.Callback = @(src,evt) sliderMoved(src,evt,hplot,hplot2,hcurPoint,voltageData,IC_dose_mean_dB,voltList); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function pulseByPulse(VolData,winNum,signalStartPoint,winInterval)
awgFreq		= VolData(1).frequency/1000;
vol			= VolData(1).voltage;
pulseCount 	= size(VolData,2);
    
subplot(4,2,6)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% FFT of each pulse at this voltage,windown is full sig range
hold on
	for i = 1:pulseCount
		sigFFT		= mag2db(VolData(i).P1);
		noiseFFT	= mag2db(VolData(i).P1noise);
		plot([VolData(1).f]',sigFFT)
		[peaks lct]	= findpeaks(sigFFT,VolData(1).f,'MinPeakProminence',35);%li's
		peaks 		= round(peaks,1);
		temp		= [reshape(lct,[size(lct,2),1]),round(reshape(lct,[size(lct,2),1])/awgFreq,1), peaks ];%li's % awgFreq awgFreq/awgFreqkhz peakValue
		temp 		= peaks;  %for print value
		plot(lct,peaks, 'm*') %labelpeaks 
		% text(lct,(peaks+8*i),num2str(round(lct'/awgFreq,1)))%li's
   end
  
plot(VolData(1).f,noiseFFT,'LineWidth',2)
xlabel('Frequency [MHz]')
ylabel('dB')
xlim([0 20])
ylim([-100 -20])
title('FFT of each pulse at this voltage,windown is full sig range')
hold off

subplot(4,2,8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% FFTs of each win in the first pulse

winInterval = winInterval/(winNum+1);
hold on 

	for i =0:winNum
		noiseWinstop   = length(VolData(1).waveforms)-1;
		noiseWinstart  = noiseWinstop - winInterval;
		signalWinStart = signalStartPoint+(winInterval*i);
		signalWinStop  = signalStartPoint+(winInterval*(i+1));
 
		pulseOrder 										= 1;% 1 the first pulse
		fourierData 									= VolData(pulseOrder).waveform(signalWinStart:signalWinStop);
		fourierDataNoise								= VolData(pulseOrder).waveform(noiseWinstart:noiseWinstop);
		[VolData(pulseOrder).P1, VolData(pulseOrder).f] = computeFFT(hanning(length(fourierData)).*fourierData, VolData(pulseOrder).sampleFreq); %
		[VolData(pulseOrder).P1noise, ~] 				= computeFFT(hanning(length(fourierDataNoise)).*fourierDataNoise, VolData(pulseOrder).sampleFreq); %
		sigFFT											= mag2db(VolData(pulseOrder).P1);
		noiseFFT									    = mag2db(VolData(pulseOrder).P1noise);
	
		plot([VolData(1).f]',sigFFT)
	
		[peaks lct]= findpeaks(sigFFT,VolData(1).f,'MinPeakProminence',15);%li's
		peaks = round(peaks,1);
		temp=[reshape(lct,[size(lct,2),1]),round(reshape(lct,[size(lct,2),1])/awgFreq,1), peaks ];%li's % awgFreq awgFreq/awgFreqkhz peakValue
		temp = peaks;  %for print value
	
		plot(lct,peaks, 'm*') %labelpeaks
		text(lct,(peaks+i*7),num2str(round(lct'/awgFreq,1)))%li's
	
	end
	
plot(VolData(1).f,noiseFFT,'LineWidth',2)
xlabel('Frequency [MHz]')
ylabel('dB')
xlim([0 20])
ylim([-100 -20])
title('FFTs of each win in the first pulse')
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end






function runAnalysis(VolData,winNum,signalStartPoint,winInterval,flag)

awgFreq			= VolData(1).frequency/1000;
vol				= VolData(1).voltage;
pulseCount 		= size(VolData,2);

noiseWinstop  	= length(VolData(1).waveforms)-1;
noiseWinstart 	= noiseWinstop-winInterval;
signalWinStart	= signalStartPoint+winInterval*winNum;
signalWinStop 	= signalStartPoint+winInterval*(winNum+1);


subplot(4,2,2)
hold on
plot(mean([VolData(1:pulseCount).waveforms],2))
ylim([-0.8 0.8])
hold on
xline(noiseWinstart,'LineWidth',2)
xline(noiseWinstop,'LineWidth',2)
xline(signalWinStart,'g','LineWidth',2)
xline(signalWinStop,'g','LineWidth',2)
title([num2str(vol) 'v']);    
        
    for i= 1:length(VolData) 
		fourierData 	 				= VolData(i).waveform(signalWinStart:signalWinStop);
		fourierDataNoise 				= VolData(i).waveform(noiseWinstart:noiseWinstop);  
		[VolData(i).P1, VolData(i).f] 	= computeFFT(hanning(length(fourierData)).*fourierData, VolData(i).sampleFreq); %
		[VolData(i).P1noise, ~] 		= computeFFT(hanning(length(fourierDataNoise)).*fourierDataNoise, VolData(i).sampleFreq); %
    end 
    
sigFFT		= mag2db(mean([VolData.P1],2)); 
noiseFFT	= mag2db(mean([VolData.P1noise],2)); 


subplot(4,2,4)
hold on
plot([VolData(1).f]',sigFFT)

[peaks lct]	= findpeaks(sigFFT,VolData(1).f,'MinPeakProminence',15);%li's
peaks 		= round(peaks,1);
temp		= [reshape(lct,[size(lct,2),1]),round(reshape(lct,[size(lct,2),1])/awgFreq,1), peaks ];%li's % awgFreq awgFreq/awgFreqkhz peakValue
temp 		= peaks;  %for print value

hold on
plot(lct,peaks, 'm*') %labelpeaks
% text(lct,peaks,[num2str(round(lct'/awgFreq,1)) num2str(peaks)])%li's

text(lct,(peaks+7*flag),num2str(round(lct'/awgFreq,1)))%li's
% text(lct+1.2,peaks, num2str(peaks) )%li's

ylim([-100 -20])
hold on
plot(VolData(1).f, noiseFFT,'LineWidth',2)
% xline(VolData(1).f(fftWinStart),'LineWidth',2)
% xline(VolData(1).f((fftWinStop)),'LineWidth',2)
xlabel('Frequency [MHz]')
ylabel('dB')
xlim([0 20])
% title('FFTs of different win')
end



function sliderMoved(src,event,hplot,hplot2,hcurPoint,voltageData,IC_dose_mean_dB,voltList)
    n = get(src,'Value');
    set(hcurPoint,'ydata',IC_dose_mean_dB(ceil(n/sum([voltageData.voltage] == voltageData(1).voltage))),'xdata',voltList(ceil(n/sum([voltageData.voltage] == voltageData(1).voltage))));
    title(num2str(round(n)))
    set(hplot,'ydata',voltageData(round(n)).waveform);
%     subplot(3,2,[2,4])
%   subplot(3,2,6)
    set(hplot2,'ydata',mag2db(voltageData(round(n)).P1));
    drawnow;

end

function [P1, f] = computeFFT(sampleData, sampleFreq)
  
L 			= (2)*length(sampleData); % li's for zero padding
rawFFT 		= fft(sampleData,L);
P2 			= abs(rawFFT/L);
P1 			= P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f 			= sampleFreq*(0:(L/2))/L;
end

function [press_y,press_y800] = DataPreProcess(voltageData) 

[~,index] 	= sortrows([voltageData.y].'); 
voltageData = voltageData(index); clear index
ts 			= voltageData;
freq  		= ts(1).frequency;
fs 			= (ts(1).sampleFreq)*1000;%li's
clear wave mean_waveform envelope startPoint endPoint press_y800

% % For 1pusle exp:  gets the pressure at each positions under same voltage but 
	for i = 1:size(voltageData,2)
		mean_waveform(i,:)	= filt(ts(i).waveforms,freq ,fs);% eliminate noise of each pulse
		envelope(i,:)   	= abs(hilbert((mean_waveform(i,:)))); %hilbert transform 
		startPoint 			= find(envelope(i,:)>max(envelope(i,:))/1.5,1,'first'); %index of first half max amplitude
		endPoint   			= find(envelope(i,:)>max(envelope(i,:))/1.5,1,'last'); %index of last half max amplitude
		Pressure(i,1)		= 1000* mean(envelope(i,(startPoint:endPoint))); %peak negative pressure
	end
	
startIndex_y800	= find([ts.y]==800,1,'first');
endIndex_y800	= find([ts.y]==800,1,'last');
press_y800 		= Pressure(startIndex_y800:endIndex_y800,1);%only output the pressure of expected position for prediction later
YList 			= unique([ts.y]);

	for j = 1:length(YList)
		startIndex	= find([ts.y]==YList(j),1,'first');
		endIndex	= find([ts.y]==YList(j),1,'last');
		press_y(j,1)= mean(Pressure((startIndex:endIndex),1));
	end


%For several pusles in 1 exp: can't get the error bar curve	 
% YList= unique([ts.y]);
% 
% 	for j = 1:length(YList)
% 		startIndex			= find([ts.y]==YList(j),1,'first');
% 		endIndex			= find([ts.y]==YList(j),1,'last');
% 		mean_waveform 		= mean(ts(startIndex:endIndex).waveforms,2);% eliminate nois
% 		envelope 			= abs(hilbert((mean_waveform))); %hilbert transform 
% 		startPoint 			= find(envelope>max(envelope)/1.5,1,'first'); %index of first half max amplitude
% 		endPoint   			= find(envelope>max(envelope)/1.5,1,'last'); %index of last half max amplitude
% 		Pressure(j,1) 		= 1000*mean(envelope(startPoint:endPoint)); %peak negative pressure
% 	end
% 	
% startIndex_y800	= find(YList==800,1,'first');
% pulseCount_y800 = find([ts.y]==800,1,'last')-find([ts.y]==800,1,'first');
% press_y800 		= Pressure(startIndex_y800,1).*ones(length(pulseCount_y800),1);
% press_y			= Pressure;

end

