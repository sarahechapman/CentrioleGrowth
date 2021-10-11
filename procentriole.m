syms kp1 kp2 kp3 kp4 kp5
syms kst1 kst2 ksas1 ksas2 kcw
syms P4(t) P4cyto(t) STIL(t) P4p(t) STILcyto(t) S6(t) S6cw(t) S6cyto(t) xa

%Defining the parameters
%These numbers are from: https://www.biorxiv.org/content/10.1101/591834v1.full
kp1 = 0.02;
kp2 = 0.02;
kp3 = 0.2;
kp4 = 0.02;
kp5 = 0.2;
kst1 = 0.16;
kst2 = 0.08;
ksas1 = 0.01;
ksas2 = 0.01;
kcw = 2;

%Defining the equations

eqn1 = diff(P4, t) ==  kp1*(1-0.1*(1/(1+exp(-0.02*(t-200)))))*(1-P4/10) - (kp2*P4+kp3*P4*(0.02+STIL^2)+kp4*P4*P4p);
eqn2 = diff(P4p, t) == (kp3*P4*(0.02+STIL^2)+kp4*P4*P4p)-kp5*P4p*(1-STIL/10);
eqn3 = diff(STIL, t) == kst1*(1/(1+exp(-0.02*(t-200))))*(0.4*P4+P4p)*(1-STIL)-kst2*STIL*(0.06+P4p)/(0.01+S6+S6cw);
eqn4 = diff(S6, t) == ksas1*(1/(1+exp(-0.1*(t-1))))*(1-S6/0.05)-ksas2*S6-kcw*P4p*STIL*S6*(1-S6cw/0.5);
eqn5 = diff(S6cw, t) == kcw*P4p*STIL*S6*(1-S6cw/0.5); 

[V] = odeToVectorField(eqn1, eqn2, eqn3, eqn4, eqn5)

M = matlabFunction(V,'vars', {'t','Y'});
sol = ode45(M,[0 1000],[0 0, 0, 0, 0]);

P4total = @(x) deval(sol,x,1)+deval(sol,x,2);
S6total = @(x) deval(sol,x,4)+deval(sol,x,5);

%Plotting the graph
fplot(P4total, [0, 1000])
hold on
fplot(S6total, [0, 1000])
fplot(@(x) deval(sol,x,3), [0,1000])
hold off
legend({'P4total','S6total','STIL'},'FontSize',20)
xlabel('Quantity (A.U.)', 'FontSize',20) 
ylabel('Time (A.U)', 'FontSize',20) 
