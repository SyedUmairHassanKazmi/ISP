clearvars -except Data
clc
Latitude = Data(1,1);
Longitude = Data(1,2);
Months = [1:12];
NOCT = 45;
HgL = Data(2,[3:14]);
TaL = Data(1,[3:14]);
Declanation_list = [-20.9, -13.0, -2.4, 9.4, 18.8, 23.1, 21.2, 13.5, 2.2, -9.6, -18.9, -23];
Htiltf = [];
Tcellf = [];
Pdcf = [];
TLMIN = -10;
TLMAX = 80;
ALMIN = -60;
ALMAX = 60;
for Loop = 1:12
    Ta = Data(1,(Loop+2));
    Hg = Data(2,(Loop+2));
    Hd = Data(3,(Loop+2));
    pg = Data(4,(Loop+2));
    Declanation = Declanation_list(Loop);
 for Tilt = TLMIN:TLMAX
  for Azimuth = ALMIN:ALMAX
   ws1 = acosd(-tand(Latitude)*tand(Declanation));
   ws2 = acosd(-tand(Latitude - Tilt)*tand(Declanation));
   ws = min([ws1 , ws2]);
   A = (cosd(Tilt) + (tand(Latitude) * cosd(Azimuth) * sind(Tilt)));
   B = (cosd(ws) * cosd(Tilt)) + (tand(Declanation) * sind(Tilt) * cosd(Azimuth));
   C = (sind(Tilt) * sind(Azimuth)) / cosd(Latitude);
   wsr = abs(min([ws , acosd(((A*B) + (C*sqrt(A^2 - B^2 + C^2)))/ (A^2 + C^2))]));
   if (A > 0 && B > 0) || (A >= B)
       wsr = -wsr;
   end
   wss = abs(min([ws , acosd(((A*B) - (C*sqrt(A^2 - B^2 + C^2)))/ (A^2 + C^2))]));
   if not((A > 0 & B > 0) | A >= B)
       wss = -wss;
   end
   a = 0.409 + 0.5016*sind(ws - 60);
   b = 0.6609 - 0.4767*sind(ws - 60);
   d = sind(ws) - (pi*ws*cosd(ws))/180;
   aa = a - (Hd/Hg);    

   if wss >= wsr
       w1 = wss;
       w2 = wsr;
       g = (1/(2*d)) * [(((0.5*b*A) - (aa * B)) * (w1 - w2) * (pi/180)) + ((((aa * A)-(b*B))*(sind(w1)-sind(w2)))-(aa*C*(cosd(w1) - cosd(w2)))) + ((0.5*b*A)*((sind(w1)*cosd(w1))-(sind(w2)*cosd(w2)))) + ((0.5*b*C)*((sind(w1)^2 - sind(w2)^2)))];
   elseif wsr > wss
       w1 = wss;
       w2 = -ws;
       g1 = (1/(2*d)) * [(((0.5*b*A) - aa * B) * (w1 - w2) * (pi/180)) + ((((aa * A)-(b*B))*(sind(w1)-sind(w2)))-(aa*C*(cosd(w1) - cosd(w2)))) + ((0.5*b*A)*((sind(w1)*cosd(w1))-(sind(w2)*cosd(w2)))) + ((0.5*b*C)*((sind(w1)^2 - sind(w2)^2)))];
       w3 = ws;
       w4 = wsr;
       g2 = (1/(2*d)) * [(((0.5*b*A) - aa * B) * (w3 - w4) * (pi/180)) + ((((aa * A)-(b*B))*(sind(w3)-sind(w4)))-(aa*C*(cosd(w3) - cosd(w4)))) + ((0.5*b*A)*((sind(w3)*cosd(w3))-(sind(w4)*cosd(w4)))) + ((0.5*b*C)*((sind(w3)^2 - sind(w4)^2)))];
       g = g1 + g2;
   end
   D = max([0 , g]);
   R = D + ((Hd/Hg) * (.5*(1+cosd(Tilt)))) + (pg*(.5*(1-cosd(Tilt))));
   Ht = (Hg * R);
   Htilt((Tilt+1+(-TLMIN)),(Azimuth+(1+(-ALMIN)))) = Ht;
   Alpha((Tilt+1+(-TLMIN)),(Azimuth+(1+(-ALMIN)))) = Azimuth;
   Beta((Tilt+1+(-TLMIN)),(Azimuth+(1+(-ALMIN)))) = Tilt;
   Tcell((Tilt+1+(-TLMIN)),(Azimuth+(1+(-ALMIN)))) = Ta + [(NOCT - 20)/800] * Ht;
   Pdc((Tilt+1+(-TLMIN)),(Azimuth+(1+(-ALMIN)))) = 0.16 * Ht * (1+((-.4/100)*(Tcell((Tilt+1+(-TLMIN)),(Azimuth+(1+(-ALMIN))))-25)));
  end
 end
Htiltf = [Htiltf Htilt];
Htiltf(imag(Htiltf) ~= 0) = 0;
Tcellf = [Tcellf Tcell];
Tcellf(imag(Tcellf) ~= 0) = 0;
Pdcf = [Pdcf Pdc];
Pdcf(imag(Pdcf) ~= 0) = 0;
end
Z = (1+(-ALMIN)+ ALMAX);
H_Jan = Htiltf(1:(1+(-TLMIN)+ TLMAX) , 1: Z);     HTo(1)  = max(H_Jan(:)); [r1 , c1] = find(H_Jan == HTo(1),  1, 'first');  To(1) = r1 - (1+(-TLMIN));   Ao(1) = c1 - (1+(-ALMIN));
H_Feb = Htiltf(1:(1+(-TLMIN)+ TLMAX), (Z+1):2*Z);    HTo(2)  = max(H_Feb(:)); [r2 , c2] = find(H_Feb == HTo(2),  1, 'first');   To(2) = r2 - (1+(-TLMIN));   Ao(2) = c2 - (1+(-ALMIN));
H_Mar = Htiltf(1:(1+(-TLMIN)+ TLMAX) , (2*Z+1):3*Z);   HTo(3)  = max(H_Mar(:)); [r3 , c3] = find(H_Mar == HTo(3),  1, 'first');  To(3) = r3 - (1+(-TLMIN));   Ao(3) = c3 - (1+(-ALMIN));
H_Apr = Htiltf(1:(1+(-TLMIN)+ TLMAX) , (3*Z+1):4*Z);   HTo(4)  = max(H_Apr(:)); [r4 , c4] = find(H_Apr == HTo(4),  1, 'first');  To(4) = r4 - (1+(-TLMIN));   Ao(4) = c4 - (1+(-ALMIN));
H_May = Htiltf(1:(1+(-TLMIN)+ TLMAX) , (4*Z+1):5*Z);   HTo(5)  = max(H_May(:)); [r5 , c5] = find(H_May == HTo(5),  1, 'first');  To(5) = r5 - (1+(-TLMIN));   Ao(5) = c5 - (1+(-ALMIN));
H_Jun = Htiltf(1:(1+(-TLMIN)+ TLMAX) , (5*Z+1):6*Z);   HTo(6)  = max(H_Jun(:)); [r6 , c6] = find(H_Jun == HTo(6),  1, 'first');  To(6) = r6 - (1+(-TLMIN));   Ao(6) = c6 - (1+(-ALMIN));
H_Jul = Htiltf(1:(1+(-TLMIN)+ TLMAX), (6*Z+1):7*Z);   HTo(7)  = max(H_Jul(:)); [r7 , c7] = find(H_Jul == HTo(7),  1, 'first');  To(7) = r7 - (1+(-TLMIN));   Ao(7) = c7 - (1+(-ALMIN));
H_Aug = Htiltf(1:(1+(-TLMIN)+ TLMAX) , (7*Z+1):8*Z);   HTo(8)  = max(H_Aug(:)); [r8 , c8] = find(H_Aug == HTo(8),  1, 'first');  To(8) = r8 - (1+(-TLMIN));   Ao(8) = c8 - (1+(-ALMIN));
H_Sep = Htiltf(1:(1+(-TLMIN)+ TLMAX) , (8*Z+1):9*Z);  HTo(9)  = max(H_Sep(:)); [r9 , c9] = find(H_Sep == HTo(9),  1, 'first');  To(9) = r9 - (1+(-TLMIN));   Ao(9) = c9 - (1+(-ALMIN));
H_Oct = Htiltf(1:(1+(-TLMIN)+ TLMAX) , (9*Z+1):10*Z); HTo(10) = max(H_Oct(:)); [r10, c10] = find(H_Oct == HTo(10),1, 'first'); To(10) = r10 - (1+(-TLMIN)); Ao(10) = c10 - (1+(-ALMIN));
H_Nov = Htiltf(1:(1+(-TLMIN)+ TLMAX) , (10*Z+1):11*Z); HTo(11) = max(H_Nov(:)); [r11, c11] = find(H_Nov == HTo(11),1, 'first'); To(11) = r11 - (1+(-TLMIN)); Ao(11) = c11 - (1+(-ALMIN));
H_Dec = Htiltf(1:(1+(-TLMIN)+ TLMAX) , (11*Z+1):12*Z); HTo(12) = max(H_Dec(:)); [r12, c12] = find(H_Dec == HTo(12),1, 'first'); To(12) = r12 - (1+(-TLMIN)); Ao(12) = c12 - (1+(-ALMIN));

for Aoo = 1:12
    if Ao(Aoo) == -60
        Ao(Aoo) = 0;
    end
end 

P_Jan = Pdcf(1:(1+(-TLMIN)+ TLMAX) , 1:Z) * 31;     Po(1)  = max(P_Jan(:)); [rp1 , cp1] = find(P_Jan == Po(1),  1, 'first');  Top(1) = rp1 - (1+(-TLMIN));   Aop(1) = cp1 - (1+(-ALMIN));                
P_Feb = Pdcf(1:(1+(-TLMIN)+ TLMAX) , (Z+1):2*Z) * 28;   Po(2)  = max(P_Feb(:)); [rp2 , cp2] = find(P_Feb == Po(2),  1, 'first');   Top(2) = rp2 - (1+(-TLMIN));  Aop(2) = cp2 - (1+(-ALMIN));
P_Mar = Pdcf(1:(1+(-TLMIN)+ TLMAX) , (2*Z+1):3*Z) * 31;   Po(3)  = max(P_Mar(:)); [rp3 , cp3] = find(P_Mar == Po(3),  1, 'first');  Top(3) = rp3 - (1+(-TLMIN));   Aop(3) = cp3 - (1+(-ALMIN));
P_Apr = Pdcf(1:(1+(-TLMIN)+ TLMAX) , (3*Z+1):4*Z) * 30;   Po(4)  = max(P_Apr(:)); [rp4 , cp4] = find(P_Apr == Po(4),  1, 'first');  Top(4) = rp4 - (1+(-TLMIN));   Aop(4) = cp4 - (1+(-ALMIN));
P_May = Pdcf(1:(1+(-TLMIN)+ TLMAX) , (4*Z+1):5*Z) * 31;   Po(5)  = max(P_May(:)); [rp5 , cp5] = find(P_May == Po(5),  1, 'first');  Top(5) = rp5 - (1+(-TLMIN));   Aop(5) = cp5 - (1+(-ALMIN));
P_Jun = Pdcf(1:(1+(-TLMIN)+ TLMAX) , (5*Z+1):6*Z) * 30;   Po(6)  = max(P_Jun(:)); [rp6 , cp6] = find(P_Jun == Po(6),  1, 'first');  Top(6) = rp6 - (1+(-TLMIN));   Aop(6) = cp6 - (1+(-ALMIN));
P_Jul = Pdcf(1:(1+(-TLMIN)+ TLMAX) , (6*Z+1):7*Z) * 31;   Po(7)  = max(P_Jul(:)); [rp7 , cp7] = find(P_Jul == Po(7),  1, 'first');  Top(7) = rp7 - (1+(-TLMIN));   Aop(7) = cp7 - (1+(-ALMIN));
P_Aug = Pdcf(1:(1+(-TLMIN)+ TLMAX) , (7*Z+1):8*Z) * 31;   Po(8)  = max(P_Aug(:)); [rp8 , cp8] = find(P_Aug == Po(8),  1, 'first');  Top(8) = rp8 - (1+(-TLMIN));   Aop(8) = cp8 - (1+(-ALMIN));
P_Sep = Pdcf(1:(1+(-TLMIN)+ TLMAX) , (8*Z+1):9*Z) * 30;  Po(9)  = max(P_Sep(:)); [rp9 , cp9] = find(P_Sep == Po(9),  1, 'first');  Top(9) = rp9 - (1+(-TLMIN));   Aop(9) = cp9 - (1+(-ALMIN));
P_Oct = Pdcf(1:(1+(-TLMIN)+ TLMAX) , (9*Z+1):10*Z) * 31; Po(10) = max(P_Oct(:)); [rp10, cp10] = find(P_Oct == Po(10),1, 'first'); Top(10) = rp10 - (1+(-TLMIN));  Aop(10) = cp10 - (1+(-ALMIN));
P_Nov = Pdcf(1:(1+(-TLMIN)+ TLMAX) , (10*Z+1):11*Z) * 30; Po(11) = max(P_Nov(:)); [rp11, cp11] = find(P_Nov == Po(11),1, 'first'); Top(11) = rp11 - (1+(-TLMIN));  Aop(11) = cp11 - (1+(-ALMIN));
P_Dec = Pdcf(1:(1+(-TLMIN)+ TLMAX) , (11*Z+1):12*Z) * 31; Po(12) = max(P_Dec(:)); [rp12, cp12] = find(P_Dec == Po(12),1, 'first'); Top(12) = rp12 - (1+(-TLMIN));  Aop(12) = cp12 - (1+(-ALMIN));

for Aoop = 1:12
    if Aop(Aoop) == -60
        Aop(Aoop) = 0;
    end
end 


Max_Decp = max(P_Dec(:));
Perc_Decp = ((Max_Decp - P_Dec)/Max_Decp)*100;
for z = 1:(1+(-TLMIN)+ TLMAX)
    for m = 1:Z
        if Perc_Decp(z,m) > 1
            Perc_Decp(z,m) = 0;
        end
    end
end
Max_Novp = max(P_Nov(:));
Perc_Novp = ((Max_Novp - P_Nov)/Max_Novp)*100;
for z = 1:(1+(-TLMIN)+ TLMAX)
    for m = 1:Z
        if Perc_Novp(z,m) > 1
            Perc_Novp(z,m) = 0;
        end
    end
end
Max_Octp = max(P_Oct(:));
Perc_Octp = ((Max_Octp - P_Oct)/Max_Octp)*100;
for z = 1:(1+(-TLMIN)+ TLMAX)
    for m = 1:Z
        if Perc_Octp(z,m) > 1
            Perc_Octp(z,m) = 0;
        end
    end
end
Max_Sepp = max(P_Sep(:));
Perc_Sepp = ((Max_Sepp - P_Sep)/Max_Sepp)*100;
for z = 1:(1+(-TLMIN)+ TLMAX)
    for m = 1:Z
        if Perc_Sepp(z,m) > 1
            Perc_Sepp(z,m) = 0;
        end
    end
end
Max_Augp = max(P_Aug(:));
Perc_Augp = ((Max_Augp - P_Aug)/Max_Augp)*100;
for z = 1:(1+(-TLMIN)+ TLMAX)
    for m = 1:Z
        if Perc_Augp(z,m) > 1
            Perc_Augp(z,m) = 0;
        end
    end
end
Max_Julp = max(P_Jul(:));
Perc_Julp = ((Max_Julp - P_Jul)/Max_Julp)*100;
for z = 1:(1+(-TLMIN)+ TLMAX)
    for m = 1:Z
        if Perc_Julp(z,m) > 1
            Perc_Julp(z,m) = 0;
        end
    end
end
Max_Junp = max(P_Jun(:));
Perc_Junp = ((Max_Junp - P_Jun)/Max_Junp)*100;
for z = 1:(1+(-TLMIN)+ TLMAX)
    for m = 1:Z
        if Perc_Junp(z,m) > 1
            Perc_Junp(z,m) = 0;
        end
    end
end
Max_Mayp = max(P_May(:));
Perc_Mayp = ((Max_Mayp - P_May)/Max_Mayp)*100;
for z = 1:(1+(-TLMIN)+ TLMAX)
    for m = 1:Z
        if Perc_Mayp(z,m) > 1
            Perc_Mayp(z,m) = 0;
        end
    end
end
Max_Aprp = max(P_Apr(:));
Perc_Aprp = ((Max_Aprp - P_Apr)/Max_Aprp)*100;
for z = 1:(1+(-TLMIN)+ TLMAX)
    for m = 1:Z
        if Perc_Aprp(z,m) > 1
            Perc_Aprp(z,m) = 0;
        end
    end
end
Max_Marp = max(P_Mar(:));
Perc_Marp = ((Max_Marp - P_Mar)/Max_Marp)*100;
for z = 1:(1+(-TLMIN)+ TLMAX)
    for m = 1:Z
        if Perc_Marp(z,m) > 1
            Perc_Marp(z,m) = 0;
        end
    end
end
Max_Febp = max(P_Feb(:));
Perc_Febp = ((Max_Febp - P_Feb)/Max_Febp)*100;
for z = 1:(1+(-TLMIN)+ TLMAX)
    for m = 1:Z
        if Perc_Febp(z,m) > 1
            Perc_Febp(z,m) = 0;
        end
    end
end
Max_Janp = max(P_Jan(:));
Perc_Janp = ((Max_Janp - P_Jan)/Max_Janp)*100;
for z = 1:(1+(-TLMIN)+ TLMAX)
    for m = 1:Z
        if Perc_Janp(z,m) > 1
            Perc_Janp(z,m) = 0;
        end
    end
end

[RT, CA] = find(Perc_Decp,1);
A_Min_Dec = CA - (1+(-ALMIN));
[RT, CA] = find(Perc_Decp,1, 'last');
A_Max_Dec = CA - (1+(-ALMIN));

[RT, CA] = find(Perc_Novp,1);
A_Min_Nov = CA - (1+(-ALMIN));
[RT, CA] = find(Perc_Novp,1, 'last');
A_Max_Nov = CA - (1+(-ALMIN));

[RT, CA] = find(Perc_Octp,1);
A_Min_Oct = CA - (1+(-ALMIN));
[RT, CA] = find(Perc_Octp,1, 'last');
A_Max_Oct = CA - (1+(-ALMIN));

[RT, CA] = find(Perc_Sepp,1);
A_Min_Sep = CA - (1+(-ALMIN));
[RT, CA] = find(Perc_Sepp,1, 'last');
A_Max_Sep = CA - (1+(-ALMIN));

[RT, CA] = find(Perc_Augp,1);
A_Min_Aug = CA - (1+(-ALMIN));
[RT, CA] = find(Perc_Augp,1, 'last');
A_Max_Aug = CA - (1+(-ALMIN));

[RT, CA] = find(Perc_Julp,1);
A_Min_Jul = CA - (1+(-ALMIN));
[RT, CA] = find(Perc_Julp,1, 'last');
A_Max_Jul = CA - (1+(-ALMIN));

[RT, CA] = find(Perc_Junp,1);
A_Min_Jun = CA - (1+(-ALMIN));
[RT, CA] = find(Perc_Junp,1, 'last');
A_Max_Jun = CA - (1+(-ALMIN));

[RT, CA] = find(Perc_Mayp,1);
A_Min_May = CA - (1+(-ALMIN));
[RT, CA] = find(Perc_Mayp,1, 'last');
A_Max_May = CA - (1+(-ALMIN));

[RT, CA] = find(Perc_Aprp,1);
A_Min_Apr = CA - (1+(-ALMIN));
[RT, CA] = find(Perc_Aprp,1, 'last');
A_Max_Apr = CA - (1+(-ALMIN));

[RT, CA] = find(Perc_Marp,1);
A_Min_Mar = CA - (1+(-ALMIN));
[RT, CA] = find(Perc_Marp,1, 'last');
A_Max_Mar = CA - (1+(-ALMIN));

[RT, CA] = find(Perc_Febp,1);
A_Min_Feb = CA - (1+(-ALMIN));
[RT, CA] = find(Perc_Febp,1, 'last');
A_Max_Feb = CA - (1+(-ALMIN));

[RT, CA] = find(Perc_Janp,1);
A_Min_Jan = CA - (1+(-ALMIN));
[RT, CA] = find(Perc_Janp,1, 'last');
A_Max_Jan = CA - (1+(-ALMIN));

A_Min_Max = [A_Min_Jan , A_Min_Feb , A_Min_Mar , A_Min_Apr , A_Min_May , A_Min_Jun , A_Min_Jul , A_Min_Aug , A_Min_Sep , A_Min_Oct , A_Min_Nov , A_Min_Dec, A_Max_Jan , A_Max_Feb , A_Max_Mar , A_Max_Apr , A_Max_May , A_Max_Jun , A_Max_Jul , A_Max_Aug , A_Max_Sep , A_Max_Oct , A_Max_Nov , A_Max_Dec];
Months2 = [1:12 , 1:12];

Mean_Tilt = mean(To);
Mean_Azimuth = mean(Ao);
Tilt_Angles = To;
Mean_Power = mean(Po);
MAX_ANNUAL_POWER = sum(Po);

figure(1);
f1 = contour(Beta,Alpha,P_Dec);
xlabel('Tilt Angles (Degrees)')
ylabel('Azimuth Angles (Degrees)')
title('December Power Contour')
hold on
plot(To(12),Ao(12), 'xm');
contour(Beta,Alpha,Perc_Decp ,'c')
legend('Power (KW-hr / m^2 / Month)', 'Max Power', '< 1% Power Difference' , 'Location', 'North' );

figure(2);
f2 = contour(Beta,Alpha,P_Jun);
xlabel('Tilt Angles (Degrees)')
ylabel('Azimuth Angles (Degrees)')
title('June Power Contour')
hold on
plot(To(6),Ao(6), 'xm');
contour(Beta,Alpha,Perc_Junp ,'c')
legend('Power (KW-hr / m^2 / Month)', 'Max Power', '< 1% Power Difference' , 'Location', 'North' );

figure(3);
f3 = plot(Months,Ao, 'r:x');
legend( f3, 'Azimuth Angles vs. Months', 'Location', 'North' );
% Label axes
xlabel('Months')
ylabel('Azimuth Angles (Degrees)')
xlim([1 12])
title('Monthly Optimum Azimuth Angle')
grid on
hold on 
plot(Months,Mean_Azimuth*ones(size(Months)), 'LineStyle','-' , 'DisplayName','Fixed Optimal Azimuth')

figure(4);

scatter(A_Min_Max , Months2 , '*r')
ylabel('Months')
xlabel('Azimuth Angles (Degrees)')
title('Min Max Azimuth Angle Values for 1% Error')
legend('Max and Min Azimuth', 'Location', 'best' );

%Figure 5
BarTemp(Months, TaL, HgL);

%Figure 6
Optimal_Tilt(Months, Tilt_Angles, Mean_Tilt);

%Figure 7
Optimal_Power(Months, Po, Mean_Power);

Final_Data = [Latitude, Longitude, Mean_Tilt, Mean_Azimuth, Mean_Power, MAX_ANNUAL_POWER, A_Min_Jun , A_Max_Jun , A_Min_Dec , A_Max_Dec ];
