%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = Example02_Domain(Demand,Arg)
  BdBox = [0 2 0 1];
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
  % Support nodes
  LeftMiddleNode = sqrt((Node(:,1)-BdBox(1)).^2+...
      (Node(:,2)-(BdBox(3)+BdBox(4))/2).^2);
  RigthMiddleNode = sqrt((Node(:,1)-BdBox(2)).^2+...
      (Node(:,2)-(BdBox(3)+BdBox(4))/2).^2);
 [foo,LeftMiddleNode] = sort(LeftMiddleNode);
 [foo,RigthMiddleNode] = sort(RigthMiddleNode);
  FixedNodes = [LeftMiddleNode(1); RigthMiddleNode(1)];
  Supp = zeros(length(FixedNodes),3);
  Supp(:,1)=FixedNodes; Supp(:,2:end)=1;
  % Load nodes
  LeftLoadNode = sqrt((Node(:,1)-0.7*(BdBox(2)/2)).^2+...
      (Node(:,2)-(BdBox(3)+BdBox(4))/2).^2);
  RigthLoadNode = sqrt((Node(:,1)-1.3*(BdBox(2)/2)).^2+...
      (Node(:,2)-(BdBox(3)+BdBox(4))/2).^2);
%   LeftLoadNode = find(abs(Node(:,1)-0.7*(BdBox(2)/2))<eps & ...
%                        abs(Node(:,2)-BdBox(4)/2)<eps);
%   RigthLoadNode = find(abs(Node(:,1)-1.3*(BdBox(2)/2))<eps & ...
%                        abs(Node(:,2)-BdBox(4)/2)<eps);                

 [foo,LeftLoadNode] = sort(LeftLoadNode);
 [foo,RigthLoadNode] = sort(RigthLoadNode);
  Load = [LeftLoadNode(1),0,50;
          RigthLoadNode(1),0,-50];
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%