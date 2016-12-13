clear
clc
close all hidden

syms r1 r2 t real
syms x y z real
syms P Pprim real
syms theta1 theta2 theta3 real
syms c1 c2 c3 real
syms s1 s2 s3 real
syms Q3 real
syms e v psi real

c1 = cos(theta1);
c2 = cos(theta2);
c3 = cos(theta3);

s1 = sin(theta1);
s2 = sin(theta2);
s3 = sin(theta3);
        
Q3      = [ c2*c3               , -c2*s3            , s2        ; ...
            c1*s3 + c3*s1*s2    , c1*c3 - s1*s2*s3  , -c2*s1    ; ...
            s1*s3 - c1*c3*s2    , c3*s1 + c1*s2*s3  ,  c1*c2    ];
        
x = r1 * cos(t);
y = r2 * sin(t);
z = 0;
P = [x ; y ; z];

Pprim = Q3 * P;









        
        