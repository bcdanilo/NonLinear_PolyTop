%-------------------------------------------------------------------------%
% Ref: H Chi, DL Ramos, AS Ramos Jr., GH Paulino, "On structural topology %
% optimization considering material nonlinearity: Plane strain versus     %
% plane stress solutions", Advances in Engineering Softwares              %
% DOI 10.1016/j.advengsoft.2018.08.017                                    %
%-------------------------------------------------------------------------%
function [z,V,fem] = NLPolyTop(fem,opt,io)
%% Initialize variables
% Optimization parameters:
Iter=0; Tol=opt.Tol*(opt.zMax-opt.zMin); Change=2*Tol; z=opt.zIni; P=opt.P;
% Threshold projection continuation parameters:
BFreq = opt.contB(1); B = opt.contB(2); 
Binc = opt.contB(3); Bmax = opt.contB(4);
% Penalization continuation parameters:
pFreq = opt.contp(1); p = opt.contp(2); 
pinc = opt.contp(3); pmax = opt.contp(4);
% Filter radius continuation parameters:
RFreq = opt.contR(1); R = opt.contR(2); 
Rinc = opt.contR(3); Rmin = opt.contR(4);
RFlg= false;
% Analysis parameters
[E,dEdy,V,dVdy] = opt.MatIntFnc(P*z,p);
% Element ID for active elements
Eid = setdiff((1:fem.NElem)',fem.Passive);
z(fem.Passive) = 1; % Set z = 1 for passive elements
% Vector to store the OF values for plotting its evolution
OF = [];
%% Plot initial density distribution
figName=[io.ProblemId,'_ElementDensities']; figure('Name',figName);
FigHandle = InitialPlot(fem,io,V); clear figName;
%% Optimization iterations
tic; 
while (Iter<opt.MaxIter) && (Change>Tol)  
  Iter = Iter + 1;
  %Compute cost functionals and analysis sensitivities
  [f,dfdE,dfdV,fem] = ObjectiveFnc(fem,E,V);
  [g,dgdE,dgdV,fem] = ConstraintFnc(fem,E,V,opt.VolFrac); 
  %Compute design sensitivities
  dfdz = P'*(dEdy.*dfdE + dVdy.*dfdV);
  dgdz = P'*(dEdy.*dgdE + dVdy.*dgdV);
  %Update design variable associated with active elements
  [z(Eid),Change] = UpdateScheme(dfdz(Eid),g,dgdz(Eid),z(Eid),opt);
  %Update material interpolation function
  if mod(Iter,BFreq)==0; B = min(B+Binc,Bmax); end
  if mod(Iter,pFreq)==0; p = min(p+pinc,pmax); end
  %Update filter radius
  if (((Change<0.05) || (Iter>=200)) && RFlg==false), RFlg= true; end
  if (mod(Iter,RFreq)==0) && (RFlg==true)
      R = max(R-Rinc,Rmin);
      P = opt.Filter(R);
  end
  %Update analysis parameters
  [E,dEdy,V,dVdy] = opt.MatIntFnc(P*z,p);
  %Output results
  fprintf(['It: %4d\tObjective: %1.4e\tChange: %1.3f\tp: %1.3f\t'...
      'r_filter: %1.4f\t||R||/||Fext||= %1.4e\tMax(|u|): %1.4e m\n'],...
      Iter,f,Change,p,R,fem.RN,max(abs([fem.U])));
  fprintf(io.fout,...
      ['It: %4d\tObjective: %1.4e\tChange: %1.3f\tp: %1.3f\t'...
      'r_filter: %1.4f\t||R||/||Fext||= %1.4e\tMax(|u|): %1.4e m\n'],...
      Iter,f,Change,p,R,fem.RN,max(abs([fem.U])));
  % Update density plots
  set(FigHandle,'FaceColor','flat','CData',1-V); drawnow
  % Store objective function value for plotting its evolution
  OF(Iter)= f;
end
t = toc;
if t<=60, fprintf('Optimization time: %i seconds \n', round(t));
elseif t<=3600, fprintf('Optimization time: %1.1f minutes \n', t/60);
elseif t<=86400, fprintf('Optimization time: %1.1f hours \n', t/3600); 
else, fprintf('Optimization time: %1.1f days \n', t/86400); 
end
fprintf(io.fout,'Optimization time: %i seconds \n', round(t));
%% Save figures
% Final topology
savefig(gcf,strcat(io.FileName,'_FinalTopology','.fig'),'compact');
% Objective function evolution
figName=[io.ProblemId,'_ObjectiveFunction']; idfig = figure('Name',figName);
axis on, box on
plot(1:length(OF),OF,'-o','Color',io.LineColor,...
    'LineWidth',io.GrLineWidth,'MarkerEdgeColor',io.MarkerColor,...
    'MarkerFaceColor',io.MarkerColor,'MarkerSize',io.MrkSz);
xlabel('Iteration'); ylabel('Objective Function');
set(gca,'fontsize',io.AxFntSz,'TickLabelInterpreter','latex');
savefig(idfig,strcat(io.FileName,'_ObjectiveFunction','.fig'),'compact');
%% ----------------------------------------------------- OBJECTIVE FUNCTION
function [f,dfdE,dfdV,fem] = ObjectiveFnc(fem,E,V)
[U,fem] = NLFEM(fem,E); %Run nonlinear FEM routine
Vext = dot(fem.Fext,U);
Wint = sum(E.*fem.Psi);
TPE = Wint - Vext; %Total potential energy
f= -TPE; % Minus sign because OF is the maximization of TPE
dfdE = -fem.Psi;
dfdV = zeros(size(V));
%% ---------------------------------------------------- CONSTRAINT FUNCTION
function [g,dgdE,dgdV,fem] = ConstraintFnc(fem,E,V,VolFrac)
g = sum(fem.ElemArea.*V)/sum(fem.ElemArea)-VolFrac;
dgdV = fem.ElemArea/sum(fem.ElemArea);
dgdE = zeros(size(E));
%% ------------------------------------------------------- VON MISES STRESS
function [VM_Stress] = von_Mises_Stress(fem,U)
V = [1 -1/2 0; -1/2 1 0; 0 0 3]; % von Mises matrix
ElemU = U(fem.eDof);    %Element displacement vectors
ee_elem = fem.B0*ElemU; %Strains at the centroid of all elements
ee_elem = reshape(ee_elem,3,[]); %Strains at the centroid of all elements
[Cauchy_S, D0] = material_model(fem.MatModel,fem.MatParam,ee_elem);
D0 = sparse(fem.rowD,fem.colD,reshape(D0,9*fem.NElem,1));
DB = D0*fem.B0;
VM_Stress = max(sqrt(sum(Cauchy_S.*(V*Cauchy_S))),eps)'; % von Mises stress
%% --------------------------------------------- OPTIMALITY CRITERIA UPDATE
function [zNew,Change] = UpdateScheme(dfdz,g,dgdz,z0,opt)  
zMin=opt.zMin; zMax=opt.zMax;  
move=opt.OCMove*(zMax-zMin); eta=opt.OCEta;
l1=0; l2=1e6;  
while l2-l1 > 1e-4
  lmid = 0.5*(l1+l2);
  B = -(dfdz./dgdz)/lmid;
  zCnd = zMin+(z0-zMin).*B.^eta;
  zNew = max(max(min(min(zCnd,z0+move),zMax),z0-move),zMin);
  if (g+dgdz'*(zNew-z0)>0),  l1=lmid;
  else                       l2=lmid;  end
end
Change = max(abs(zNew-z0))/(zMax-zMin);
%% ----------------------------------------------------------- INITIAL PLOT
function [handle] = InitialPlot(fem,io,z01)
ElemNodes = cellfun(@length,fem.Element); %Number of nodes of each element
Faces = NaN(fem.NElem,max(ElemNodes));    %Populate Faces with NaN
for el = 1:fem.NElem; Faces(el,1:ElemNodes(el)) = fem.Element{el}(:); end
title({'Element Densities';''},'FontSize',io.TtlFntSz);
patch('Faces',Faces,'Vertices',fem.Node,'FaceVertexCData',0.*z01,...
      'FaceColor','flat','EdgeColor','k','linewidth',1.5);
handle= patch('Faces',Faces,'Vertices',fem.Node,'FaceVertexCData',...
                1-z01,'FaceColor','flat','EdgeColor','none');
axis equal; axis off; axis tight; colormap(gray); caxis([0 1]);
drawnow;
%-------------------------------------------------------------------------%
%------------------------ NLPolyTop - History ----------------------------%
%
% history: Created:    13-05-2021   Danilo Cavalcanti
%          Supervised by:           Sylvia R. M. Almeida
%
% The routines were created based on the original PolyTop and on PolyStress
% adapting them whenever was necessary
%-------------------------------------------------------------------------%