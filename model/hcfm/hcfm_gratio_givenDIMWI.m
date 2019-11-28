
function g = hcfm_gratio_givenDIMWI(fitres,imgPara)


Amw = fitres.S0_MW;
Aiw = fitres.S0_IEW .* imgPara.icvf;

g = sqrt(abs(Aiw)./(abs(Aiw)+abs(Amw)/imgPara.rho_mw));


end