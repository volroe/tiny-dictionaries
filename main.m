%% setup parameter
errorthresh=0.06;
load('theta','theta_t');
TR=5E-3;
B1slice=1;

T1=1.088;
T2=0.069;
Rho=1;

%% generate dictionary
[dict_adapt,lookup,jacobians] = createDictionary_recursion(errorthresh,theta_t,TR,B1slice);

%% create probing signal
[signals_probe,lookup_probe,jac_probe] = HSFP_signal([T1,T2],theta_t,B1slice,TR,[0 10 0 10]);
signals_probe(:,1)=[];%not measured
jac_probe(:,1,:)=[];%not measured

%% project probing signal onto manifold using generated dictionary
[param,x_pr, x_match, y_match, PD]=projectToManifold(signals_probe,lookup,dict_adapt,jacobians);

%% print out results
disp("--------------")
fprintf('ground truth   \n T1: %f T2:%f PD:%f\n\n',T1,T2,Rho)
fprintf('reconstruction \n T1: %f T2:%f PD:%f\n',param(1),param(2),param(3))