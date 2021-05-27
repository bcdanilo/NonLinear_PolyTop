%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = Example01_Domain(Demand,Arg)
  BdBox = [0 1 0 1];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  Dist = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  LeftMiddleNode = find((abs(Node(:,1)-BdBox(1))<eps) & ...
                         abs(Node(:,2)-BdBox(4)/2)<eps);
  MiddleNode = find(abs(Node(:,1)-BdBox(2)/2)<eps & ...
                    abs(Node(:,2)-BdBox(4)/2)<eps);
  RigthMiddleNode = find(abs(Node(:,1)-BdBox(2))<eps & ...
                         abs(Node(:,2)-BdBox(4)/2)<eps);
  FixedNodes = [LeftMiddleNode; RigthMiddleNode];
  Supp = zeros(length(FixedNodes),3);
  Supp(:,1)=FixedNodes; Supp(:,2:end)=1; 
  Load = [MiddleNode,0,-1000];
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%