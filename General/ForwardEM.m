function[QP] = ForwardEM(~,param)

S.x=0;                      %x-coordinate receiver (m)
S.y=0;                      %y-coordinate receiver (m)
S.z=0;                      %z-coordinate receiver (m)
S.height=0.16;              %Height of transmitter (m)
S.freq=9000;                %Frequency (Hz)
S.mom=1;                    %Transmitter moment (A.m^2)
S.r=1;                      %Coil spacing (m)
S.ori='ZZ';                 %Coil orientation (2 letter combination of X, Y, and Z)

tmp = length(param);
nb_layer = (tmp+1)/3;

M.sus=linspace(100,100,nb_layer) .* 1e-5;   %Susceptibility of layer(s) (SI unit)
M.con=param(nb_layer:2*nb_layer-1);         %Conductivity of layer(s) (S/m)
M.perm=param(nb_layer:end);                 %Permittivity of layer(s) (F/m)
M.thick=[param(1:nb_layer-1) 100];          %Layer(s) thickness (m)


[~,QP] = FDEM1DFWD_RC(S,M);

end