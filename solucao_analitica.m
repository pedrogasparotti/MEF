function y = solAnalitica(xPontos)
global E L I q0 F0 M0 cf

% Constantes auxiliares
r = (sqrt(2)/2)*((cf/(E*I))^(1/4));
A = sin(r*L) / cos(r*L);
B = -(M0*exp(r*L))/(2*(r^2)*E*I*cos(r*L));
G = (F0*exp(r*L))/(2*(r^3)*E*I*cos(r*L));

% Defini��o das Constantes para resolver o sistema da EDO
syms C1 C2 C3 C4

% Equa��es das Condi��es de Contorno Aplicadas
eqn1 = C1 + C3 == -q0/cf;
eqn2 = C1 + C2 - C3 + C4 == 0;
eqn3 = -A*C1*exp(2*r*L) + C2*exp(2*r*L) + A*C3 - C4 == B;
eqn4 = (1+A)*C1*exp(2*r*L) - (1-A)*C2*exp(2*r*L) - (1-A)*C3 ...
       - (1+A)*C4 == G;
% Vetor com as quatro equa��es
eqns = [eqn1 eqn2 eqn3 eqn4];

% Vetor com as inc�gnitas do sistema
vars = [C1 C2 C3 C4];

% Resolu��o do Sistema
[c1 ,c2, c3, c4] = solve(eqns, vars);

% Transforma a solu��o simb�lica para real
c1 = eval(c1);
c2 = eval(c2);
c3 = eval(c3);
c4 = eval(c4);

% C�lculo do valor da fun��o aproximada
for i = 1:size(xPontos,2)
    f1 = exp(r*xPontos(i))*(c1*cos(r*xPontos(i))+c2*sin(r*xPontos(i)));
    f2 = exp(-r*xPontos(i))*(c3*cos(r*xPontos(i))+c4*sin(r*xPontos(i)));
    y(i) = f1 + f2 + (q0/cf);
end
