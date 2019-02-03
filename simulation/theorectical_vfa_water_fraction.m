Amy = 0.1;
Aex = 0.3;
Aax = 0.6;

t1my = 118e-3;
t1ax = 1;
t1ex = 1;

tr = 46e-3;

fa_r = [0:90];

E1my = exp(-tr./t1my);
E1ax = exp(-tr./t1ax);
E1ex = exp(-tr./t1ex);
for kfa = 1:length(fa_r)
    fa = fa_r(kfa);

smw(kfa) = (Amy.*(1-E1my).*sind(fa))./(1-E1my.*cosd(fa));
siw(kfa) = (Aax.*(1-E1ax).*sind(fa))./(1-E1ax.*cosd(fa));
sew(kfa) = (Aex.*(1-E1ex).*sind(fa))./(1-E1ex.*cosd(fa));
end
subplot(131);plot(smw./(smw+siw+sew),'x-');title('MWF');xlabel('FA');ylabel('%');
subplot(132);plot(siw./(smw+siw+sew),'x-');title('IWF');xlabel('FA');
subplot(133);plot(sew./(smw+siw+sew),'x-');title('EWF');xlabel('FA');

% export_fig theorectical_water_fraction -png -r300

%%

