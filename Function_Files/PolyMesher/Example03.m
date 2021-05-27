%----------------- PolyTop with Material Nonlinearity --------------------%
% Author: Danilo Cavalcanti (Federal University of Goiás)                 %
% Implementation was made using PolyTop v1.1 and PolyStress v1.1          %
%                                                                         %
% Ref: H Chi, DL Ramos, AS Ramos Jr., GH Paulino, "On structural topology %
% optimization considering material nonlinearity: Plane strain versus     %
% plane stress solutions", Advances in Engineering Softwares              %
% DOI 10.1016/j.advengsoft.2018.08.017                                    %
%-------------------------------------------------------------------------%
function [Node,Element,Supp,Load] = Example03(ProblemId)
filename= ['./Examples/' ProblemId '/' ProblemId '_Data.mat'];
if exist(filename,'file')==0
    [Node,Element,Supp,Load,~] = PolyMesher(@Example03_Domain,40000,50);
    save(filename,'Node','Element','Supp','Load');
else
    load(filename);
    PolyMshr_PlotMsh(Node,Element,size(Element,1),Supp,Load);
end
function PolyMshr_PlotMsh(Node,Element,NElem,Supp,Load)
clf; axis equal; axis off; hold on;
Element = Element(1:NElem)';                 %Only plot the first block
MaxNVer = max(cellfun(@numel,Element));      %Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
ElemMat = cellfun(PadWNaN,Element,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
patch('Faces',ElemMat,'Vertices',Node,'FaceColor','w'); pause(1e-6)
idx=1;
if exist('Supp','var')&&~isempty(Supp) %Plot Supp BC if specified
  h(idx)=plot(Node(Supp(:,1),1),Node(Supp(:,1),2),'b>','MarkerSize',8);
  leg{idx}='Support';
  idx=idx+1;
end
if exist('Load','var')&&~isempty(Load) %Plot Load BC if specified
  h(idx)=plot(Node(Load(:,1),1),Node(Load(:,1),2),'m^','MarkerSize',8);
  leg{idx}='Load';
  legend(h,leg,'Location','Northoutside','Orientation','Horizontal')
end
%-------------------------------------------------------------------------%