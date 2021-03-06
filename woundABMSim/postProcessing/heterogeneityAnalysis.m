function [AveFiberDotDistances,AveFiberDotMVL] = heterogeneityAnalysis(FiberInfoStruct)
% This function performs the heterogeneity alaysis done in Will
% Richardson's 2016 paper in the Biophysical Journal.

% INPUT:
%   FiberAngle: array of vectors containing the angles of all grids in the
%   infarct for separate model runs
%   FiberPosX: array of vectors containing the x positions corresponding to
%   the above for separate model runs
%   FiberPosY: array of vectors containing the y position corresponding to
%   the above for separate model runs

% ANALYSIS
FiberAngle = FiberInfoStruct.FiberAngle;
FiberPosX = FiberInfoStruct.FiberPosX;
FiberPosY = FiberInfoStruct.FiberPosY;

% data = [2.63722189	5.545375549	-0.433640291	-7.163005522	2.904619959	6.896319885	12.82131763	11.29678547	14.01065201	66.70786968	60.95674887	13.32014657	-24.408353	-57.15619854	-9.732409206	-7.601026157	17.34884718	9.380660914	15.30725933	63.16980206	-17.65123838	44.37177687	42.07797445	8.484195597	-54.19169696	-41.12537951	5.466286295	-29.90632034	29.79933184	10.63653074	27.33175212	12.43250487	87.19210406	32.0081503	-87.83459189	60.57429562	55.11651828	89.56416627	82.83110624	64.52949388	23.67502196	39.32438819	-57.94047061	-87.31634105	6.694642976	30.81103893	38.29940782	-60.29035343	-45.7315672	-0.014532717	4.233590852	-56.39476544	-30.85970961	-43.72963916	32.24112582	28.36523785	-85.11363754	49.55538101	76.64803632	84.91521455	-55.7679536	13.69215328	32.55056361	77.65922074	-8.648087784	77.66538682	-54.5661013	88.59042281	36.43436029	12.28718166	48.38832447	-5.363000542	5.770648838	21.07051569	3.273562384	-61.60453412	66.30605144	67.67920107	-76.57917444	-24.0494653	17.05284027	43.580024	-14.81084916	-50.88120938	88.44344541	43.65482153	-23.02971124	-39.40301628	-16.48851098	-14.97887634	-4.544186943	-80.10302367	-81.61132688	-45.85997419	36.38754784	-35.85971061	-42.84375063	-24.48040837	-5.486801609	-40.07948184	-57.08917926	24.82812548	-85.94466136	-37.30671552	-74.37866653	-64.89405821	-82.31694026	-81.62965489	-23.97671815	54.86889344	11.57081773	-71.16178813	-27.25196795	45.78432452	-64.51560479	81.63990653	-72.27086164	58.44348625	24.39957441	45.89148113	17.70285768	-63.22209626	37.33056659	-65.7772111	-57.05880409	86.32344435	-7.752882161	77.63554137	-86.68382823	46.17388704	-67.61588616	28.00392988	-45.29284852	-3.985903714	44.22150054	-16.4490979	-60.55909961	21.6244016	-9.160157306	43.70086147	-67.35448679	-17.21767008	-55.23861562	13.55088905	-72.20597189	-70.73045613	0.324807686	-59.48597612	76.04925405	-40.12568932	-69.95377131	-51.59938439	63.53030354	42.90518243	-89.45216434	-59.84264771	-82.04956674	-40.96925937	-77.3293267	-0.618587476	-63.73704406	86.42303503	-40.79698938	15.22434063	-65.08756431	-59.12373352	18.92589464	74.40075873	27.08731471	82.106393	54.80163978	-46.49974737	47.11198484	21.63829036	12.77660559	-49.99723961	-55.47631572	-73.98709583	11.18278868	-48.67431368	-19.85571049	-40.60418981	2.03878694	71.86043404	49.60340525	47.70267311	-55.43981704	54.71334084	60.54897629	-5.05149734	71.43939443	-69.14742532	64.60566841	-60.87573365	78.08534533	-62.47337615	70.04334068	49.23576741	36.61574672	36.23678385	7.007744099	-78.17390978	-33.77446585	78.73193191	15.3945515	65.22958154	-71.18786437	12.09194215	63.47621893	85.59937642	-61.69314023	47.60342957	-20.42472215	-76.95353615	84.83820896	-41.70103883	-87.20319087	-19.54336942	-73.46265233	-63.5020929	17.1464341	4.320016351	-21.00645233	-29.88352773	-34.98683689	68.56485209	76.61852966	-0.506084653	60.12693241	36.3542941	36.22743065	-89.32441612	-63.65095104	79.59988926	-77.39924688	-13.14698673	68.39664041	16.10230138	-40.819092	40.42035487	54.02849879	-17.77218433	24.19412725	66.49138241	21.74793596	69.39256639	0.291692885	27.03491846	56.78144768	-79.68397096	0.6706762	72.42311424	-14.57700939	-83.84538798	-47.46775717	43.3805399	-89.66302228	48.34721917	-74.55790285	-40.99731445	73.24232598	-57.96279947	-26.66718768	62.51867688	-64.14224877	-86.79839227	-57.48209943	-49.21234324	60.26985502	-8.827849781	13.81264841	0.082641147	0.383790248	30.72318543	0.221371521	38.91195705	24.71479063	-85.06618513	30.56171649	-33.33931424	86.28352756	70.40626935	-75.97408032	-42.01351101	81.00141465	-56.61452436	42.43887646	87.26707602	-64.51074507	-45.76335395	-6.93E+01	51.83680648	-54.8791374	-22.696462	41.19800678	-30.63477337	63.29119091	12.41111262	-71.48004538	84.71468696	-39.85315708	-75.90146224	61.32709632	74.97156298	-67.92854089	1.824243811	72.32250899	-71.38375633	56.88434875	11.87241341	8.089485259	84.67557168	11.71332946	-4.222277455	11.66067627	-12.13396138	10.71022988	-7.382663559	63.25208989	25.01366888	32.65330415	-43.20174828	-72.46363344	88.54188665	-13.94375321	45.98239725	35.04362806	73.92376159	57.17972249	-65.4980504	-75.4672666	-68.03070753	-40.30664028	-18.92895321	-54.17334322	-60.3069117	63.93227342	63.85058329	-3.325327299	-3.307662754	-45.98985008	-87.76579806	-0.910772697	-17.51769555	-84.56333341	-35.02676434	-62.32209291	-2.344477017	50.37846934	5.041307468	77.59907106	-32.76346781	-5.140825095	-16.77833901	14.18137065	72.75736068	35.53144838	5.221089587	-53.13176307	9.739713753	63.62676827	-5.060028999	28.87677947	75.22400438	40.89359081	43.54925059	54.2434611	17.65888432	-35.43441352	-20.00931257	-32.74490674	-38.83772299	-63.57125649	-40.71045554	52.90783699	81.02741494	-19.71854163	56.66422165	82.38625176	40.57371635	-25.39122989	-0.370870835	-58.74625926	-68.33299509	67.5125649	-21.61784349	4.814346333	-21.37462596	-86.36519814	14.91062153	-0.631035973	-32.01413671	-28.29784122	-13.46285649	65.42936244	67.30783492	10.58360714	-33.22319441	26.19354727	27.9284264	71.79908452;
%         13	13	13	13	13	14	14	14	14	14	14	14	14	14	14	14	15	15	15	15	15	15	15	15	15	15	15	15	15	16	16	16	16	16	16	16	16	16	16	16	16	16	16	16	17	17	17	17	17	17	17	17	17	17	17	17	17	17	17	17	17	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	19	19	19	19	19	19	19	19	19	19	19	19	19	19	19	19	19	19	19	19	19	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	21	21	21	21	21	21	21	21	21	21	21	21	21	21	21	21	21	21	21	21	21	22	22	22	22	22	22	22	22	22	22	22	22	22	22	22	22	22	22	22	22	22	22	22	23	23	23	23	23	23	23	23	23	23	23	23	23	23	23	23	23	23	23	23	23	23	23	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	26	26	26	26	26	26	26	26	26	26	26	26	26	26	26	26	26	26	26	26	26	26	26	27	27	27	27	27	27	27	27	27	27	27	27	27	27	27	27	27	27	27	27	27	28	28	28	28	28	28	28	28	28	28	28	28	28	28	28	28	28	28	28	28	28	29	29	29	29	29	29	29	29	29	29	29	29	29	29	29	29	29	29	29	29	29	30	30	30	30	30	30	30	30	30	30	30	30	30	30	30	30	30	30	30	31	31	31	31	31	31	31	31	31	31	31	31	31	31	31	31	31	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	33	33	33	33	33	33	33	33	33	33	33	33	33	34	34	34	34	34	34	34	34	34	34	34	35	35	35	35	35;
%         22	23	24	25	26	19	20	21	22	23	24	25	26	27	28	29	18	19	20	21	22	23	24	25	26	27	28	29	30	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	18	19	20	21	22	23	24	25	26	27	28	29	30	19	20	21	22	23	24	25	26	27	28	29	22	23	24	25	26];
% FiberAngle = data(1,:)';
% FiberPosX = data(2,:)';
% FiberPosY = data(3,:)';
    
NumSims = size(FiberAngle,3); % number of simulation replicates for averaging results
T = 1; % vector of timepoints to be analyzed

for s = 1:NumSims
    for t = 1:1:length(T)
      
        % fibers
        fibersample = min(5000,length(FiberAngle));
        fiberdotmatave = [];
        fiberdotmatstdev = [];
        fiberdotmatave2 = [];
        fiberdotmatstdev2 = [];
        fiberdotmatmixave = [];
        fiberdotmatmixstdev = [];

        testfibers1 = randperm(length(FiberAngle),fibersample);
        testfibers2 = randperm(length(FiberAngle),fibersample);
        FiberPosX1 = repmat(FiberPosX(testfibers1),1,length(testfibers1));
        FiberPosY1 = repmat(FiberPosY(testfibers1),1,length(testfibers1));
        FiberAngle1 = repmat(FiberAngle(testfibers1,t,s),1,length(testfibers1));
        FiberPosX2 = repmat(FiberPosX(testfibers2),1,length(testfibers2));
        FiberPosY2 = repmat(FiberPosY(testfibers2),1,length(testfibers2));
        FiberAngle2 = repmat(FiberAngle(testfibers2,t,s),1,length(testfibers2));

        fiberdistances = sqrt((FiberPosX1-FiberPosX2').^2 + (FiberPosY1-FiberPosY2').^2);
        fiberdistances = reshape(fiberdistances,numel(fiberdistances),1);
        fiberCbar = (cos(2*FiberAngle1*pi/180) + cos(2*FiberAngle2'*pi/180))/2;
        fiberSbar = (sin(2*FiberAngle1*pi/180) + sin(2*FiberAngle2'*pi/180))/2;
        fiberalignment = sqrt(fiberCbar.^2 + fiberSbar.^2);
        fiberalignment = reshape(fiberalignment,numel(fiberalignment),1);
        
        FiberDotMat = sortrows([fiberdistances(fiberdistances>0) fiberalignment(fiberdistances>0)]);
        shufflei = randperm(size(FiberDotMat,1));
        FiberDotMatmixed = [FiberDotMat(:,1) FiberDotMat(shufflei,2)];

        DistBins = [0:max(fiberdistances)/100:max(fiberdistances)];
        for k=2:length(DistBins)
            rangelow = DistBins(k-1);
            rangehigh = DistBins(k);
            indexbin = find(FiberDotMat(:,1) > rangelow & FiberDotMat(:,1) <= rangehigh);
            if numel(indexbin) >1
                fiberdotmatave = [fiberdotmatave; mean(FiberDotMat(indexbin,:))];
                fiberdotmatstdev = [fiberdotmatstdev; std(FiberDotMat(indexbin,:))];
                fiberdotmatmixave = [fiberdotmatmixave; mean(FiberDotMatmixed(indexbin,:))];
                fiberdotmatmixstdev = [fiberdotmatmixstdev; std(FiberDotMatmixed(indexbin,:))];
            end
        end
        
        FiberHTRGNTYpeak(t,s) = fiberdotmatave(1,2);
        meanMVL = mean(fiberdotmatave(:,2));
        fiberheterogeneity = fiberdotmatave(:,2) - meanMVL;
        fiberintersection = min(find(fiberheterogeneity <= 0));
        if fiberintersection ==1
            FiberHTRGNTYdist(t,s) = 0;
            FiberHTRGNTYarea(t,s) = 0;
        else
            FiberHTRGNTYdist(t,s) = fiberdotmatave(fiberintersection,1);
            FiberHTRGNTYarea(t,s) = trapz(fiberdotmatave(1:fiberintersection,1),fiberheterogeneity(1:fiberintersection));
        end
        
        if s==1 && t==1
            FiberDotDistances(:,t,s) = fiberdotmatave(:,1);
            FiberDotMVL(:,t,s) = fiberdotmatave(:,2);
        elseif size(fiberdotmatave,1) > size(FiberDotDistances,1)
            FiberDotDistances(:,t,s) = fiberdotmatave(1:size(FiberDotDistances,1),1);
            FiberDotMVL(:,t,s) = fiberdotmatave(1:size(FiberDotDistances,1),2);
        else % when size(fiberdotmatave,1) <= size(FiberDotDistances,1)
            FiberDotDistances = FiberDotDistances(1:size(fiberdotmatave,1),:,:);
            FiberDotMVL = FiberDotMVL(1:size(fiberdotmatave,1),:,:);
            FiberDotDistances(:,t,s) = fiberdotmatave(:,1);
            FiberDotMVL(:,t,s) = fiberdotmatave(:,2);
        end   
    end
end

AveFiberDotDistances = mean(FiberDotDistances,3);
AveFiberDotMVL = mean(FiberDotMVL,3);
% AveFiberHTRGNTYpeak = mean(FiberHTRGNTYpeak,2);
% AveFiberHTRGNTYdist = mean(FiberHTRGNTYdist,2);
% AveFiberHTRGNTYarea = mean(FiberHTRGNTYarea,2);
% for t = 1:1:length(T)
%     meanMVL = mean(AveFiberDotMVL(:,t));
%     FiberHTRGNTYpeakAve(t) = AveFiberDotMVL(1,t);
%     fiberheterogeneity = AveFiberDotMVL(:,t) - meanMVL;
%     fiberintersection = min(find(fiberheterogeneity <= 0));   
%     if fiberintersection ==1
%         FiberHTRGNTYdistAve(t) = 0;
%         FiberHTRGNTYareaAve(t) = 0;
%     else
%         FiberHTRGNTYdistAve(t) = FiberDotDistances(fiberintersection,t);
%         FiberHTRGNTYareaAve(t) = trapz(FiberDotDistances(1:fiberintersection,t),fiberheterogeneity(1:fiberintersection));
%     end
% end
     
% StdevFiberDotMVL = std(FiberDotMVL,1,3);
% StdevFiberHTRGNTYpeak = std(FiberHTRGNTYpeak,1,2);
% StdevFiberHTRGNTYdist = std(FiberHTRGNTYdist,1,2);
% StdevFiberHTRGNTYarea = std(FiberHTRGNTYarea,1,2);

%plot(AveFiberDotDistances(:,end),AveFiberDotMVL(:,end))