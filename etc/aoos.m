function Z_av=aoos(Z,featln,N)
% this funcion will do averaging over one sample
% featln is feature length
% N is original tesing length
% Z_av is averaged per sample

M_Z=1/featln*Z*blockones(N/featln,featln);
Z_av=M_Z(:,1:featln:end);

end %end of the function file