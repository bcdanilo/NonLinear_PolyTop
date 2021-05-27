%----------------- PolyTop with Material Nonlinearity --------------------%
% Author: Danilo Cavalcanti (Federal University of Goiás)                 %
% Implementation was made using PolyTop v1.1 and PolyStress v1.1          %
%                                                                         %
% Ref: H Chi, DL Ramos, AS Ramos Jr., GH Paulino, "On structural topology %
% optimization considering material nonlinearity: Plane strain versus     %
% plane stress solutions", Advances in Engineering Softwares              %
% DOI 10.1016/j.advengsoft.2018.08.017                                    %
%-------------------------------------------------------------------------%
clear;close all;fclose all;clc; 
restoredefaultpath; addpath(genpath('./')); %Use all folders and subfolders
set(0,'defaulttextinterpreter','latex');
fprintf('**********************************************************\n');
fprintf('* Nonlinear PolyTop: Considering material nonlinearity   *\n');
fprintf('* OF: Maximizing stationary total potential energy       *\n');
fprintf('* Volume constraint and single material                  *\n');
fprintf('* Created by: Danilo Cavalcanti            May 2021      *\n');
fprintf('* Supervised by: Sylvia Almeida                          *\n');
fprintf('**********************************************************\n\n');
%% ------------------------------------------------------------ CREATE Mesh
fprintf('*** Mesh generation with PolyMesher:\n');
%Example 01: Simply supported square domain with central load -------------
% ProblemId= 'Example01'; [~,~] = mkdir(['./Examples'],ProblemId);
% [Node,Element,Supp,Load] = PolyMesher(@Example01_Domain,40000,0,P);

%Example 02: Simply supported rectangle domain with opposite loads --------
% ProblemId= 'Example02'; [~,~] = mkdir(['./Examples'],ProblemId);
% [Node,Element,Supp,Load] = PolyMesher(@Example02_Domain,40000,50);

%Example 03: Top central load in a laterally constrained rectangular domain
ProblemId= 'Example03_P5000_PlaneStress'; 
[~,~] = mkdir(['./Examples'],ProblemId);
[Node,Element,Supp,Load] = Example03(ProblemId);

%Test: Internally pressurized cylinder (quarter of the domain)
% [Node,Element,Supp,Load] = PolyMesher(@QuarterCylinder_Domain,1200,30);
% return

NElem = size(Element,1); % Number of elements

% Save figure with mesh and boundary conditions
filename= ['./Examples/' ProblemId '/' ProblemId];
savefig(gcf,strcat(filename,'_Mesh_and_BC','.fig'),'compact');

%% ---------------------------------------------------- CREATE 'fem' STRUCT
E0 = 70e6; % E0 in kPa 
nu = 0.3;
%G = E0/2.5; Et = E0; Ec = 0.1*E0;  % 0<=(Et,Ec)<=3*G; %Material props. (linear)
alpha = 1500;
beta = nu/(1-2*nu);
mu = E0/alpha/(1+nu);
fem = struct(...
  'NNode',size(Node,1),...      % Number of nodes
  'NElem',size(Element,1),...   % Number of elements
  'Node',Node,...               % [NNode x 2] array of nodes
  'Element',{Element},...       % [NElement x Var] cell array of elements
  'Supp',Supp,...               % Array of supports
  'Load',Load,...               % Array of loads
  'Passive',[],...              % Passive elements  
  'Thickness',1,...             % Element thickness
  'MatModel','Ogden',...        % Material model ('Bilinear','Ogden')
  'MatParam',[mu,alpha,beta],...      % Material parameters for MatModel
  'MatCase','PlaneStress',...
  'TolR', 1e-4, ...             % Tolerance for norm of force residual
  'MaxIter', 30, ...            % Max NR iterations per load step
  'MEX', 'Yes');                % Tag to use MEX functions in NLFEM routine
%% ---------------------------------------------------- CREATE 'opt' STRUCT
R0 = 0.04; q = 1; % Filter radius and filter exponent
VolFrac = 0.15;
p = 1;
m = @(y,p)MatIntFnc(y,'SIMP',[p]);
filter = @(R)PolyFilter(fem,R,q);
%P = PolyFilter(fem,R0,q);
P = filter(R0);
zIni = 0.5*ones(size(P,2),1);
opt = struct(...               
  'zMin',0.0,...              % Lower bound for design variables
  'zMax',1.0,...              % Upper bound for design variables
  'zIni',zIni,...             % Initial design variables
  'MatIntFnc',m,...           % Handle to material interpolation fnc.
  'contp',[2,1,0.1,3],...     % Penalization (SIMP) continuation params. 
  'contB',[2,1,0,1],...       % Threshold projection continuation params.
  'contR',[10,R0,R0/4,R0/4],...% Filter radius continuation params.
  'Filter',filter,...         % Handle to PolyFilter to use in cont. radius
  'P',P,...                   % Matrix that maps design to element vars.
  'VolFrac',VolFrac,...       % Specified volume fraction constraint
  'Tol',0.01,...              % Convergence tolerance on design vars.
  'MaxIter',300,...           % Max. number of optimization iterations
  'OCMove',0.2,...            % Allowable move step in OC update scheme
  'OCEta',0.5 ...             % Exponent used in OC update scheme
   );
%% ----------------------------------------------------- CREATE 'io' STRUCT
fout= fopen(strcat(filename,'.out'),'w');   % Open utput data file
io= struct(...
  'ProblemId',ProblemId,...     % Problem identification label
  'FileName',filename,...       % Path + Initial name of the file
  'fout',fout,...               % Output data file id
  'LineColor','r',...           % Color of the lines in the final graphs
  'MarkerColor','b',...         % Color of the markers in the final graphs
  'GrLineWidth',2,...           % Line width for the graphs
  'MrkSz',5,...                 % Marker size
  'TtlFntSz',16,...             % Title font size
  'AxFntSz',16,...              % Axes numbers font size
  'AxInt','latex',...           % Axes label interpreter ('txt' or 'latex')
  'LblFntSz',16 ...             % Axes names font size
   );
%% ---------------------------------------------------------- RUN 'PolyTop'
fem = preComputations(fem); % Run preComputations before running NLPolyTop
fprintf('\n*** Start optimization iterations:\n');
[z,V,fem] = NLPolyTop(fem,opt,io);

% for penal = 1:0.5:4        %Continuation on the penalty parameter
%    fprintf(['current p: ', num2str(penal)]);
%    opt.MatIntFnc = @(y)MatIntFnc(y,'SIMP',penal);
%    [opt.zIni,V,fem] = PolyTop(fem,opt);
% end
save([filename,'_Result.mat'],'fem','z','V');
%-------------------------------------------------------------------------%