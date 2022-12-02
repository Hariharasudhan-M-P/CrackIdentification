close all
clear
clc

crackRatio = [0,.2,.4,.6,.8];
naturalFrequencys = [];
for kkk = 1:5
  crackLength = crackRatio(kkk)*.5;
  area  = .5*.5;
  totalLength = 30;
  numberOfNodes = 30;
  rho = 7850;
  numberOfElements = numberOfNodes-1;
  moiCrack = ((.5-crackLength)^4)/12;
  elasticityModulus = 2*10^11;
  I = (.5^4)/12;
  KG = zeros(2*numberOfNodes,2*numberOfNodes);
  M = zeros(2*numberOfNodes,2*numberOfNodes);
  for i=1:numberOfElements
    icon(i,1) = i;
    icon(i,2) = i+1;
  end

  for lmn = 1:numberOfElements
    n1 = icon(lmn,1);
    n2 = icon(lmn,2);
    le = totalLength/numberOfElements;

    k = elasticityModulus*I/le^3;
    if lmn >= 10 && lmn <= 11
      k = elasticityModulus*(moiCrack)/le^3;
    endif
    i1 = icon(lmn,1); i2 = icon(lmn,2);
    kl = k * [12 6*le -12 6*le;
              6*le 4*le^2 -6*le 2*le^2;
              -12 -6*le 12 -6*le;
              6*le 2*le^2 -6*le 4*le^2];
    iv = [2*i1-1 2*i1 2*i2-1 2*i2];
    m=(rho*area*le/420)*[156, 22*le, 54, -13*le;
                      22*le, 4*le^2, 13*le, -3*le^2;
                      54, 13*le, 156, -22*le;
                     -13*le, -3*le^2, -22*le, 4*le^2];
    M(iv,iv) += m;
    KG(iv,iv) += kl;
  end
  ibcV = 1;
  ibcV2 = 2*numberOfNodes-1;
  bcV = 0;
  bcTheta = 0;
  KG(:,ibcV) = bcV;
  KG(ibcV,:) = bcTheta;
  KG(:,ibcV2) = bcV;
  KG(ibcV2,:) = bcTheta;
  KG(1,1) = 1;
  KG(ibcV2,ibcV2) = 1;
  cracked = sqrt(eig(KG,M));
  naturalFrequencys(kkk) = cracked(4);
end
crackRatio'
naturalFrequencys'