for L = 1:12
    Ta = Data(1,(L+2));
    h = Data(2,(L+2));
    hd = Data(3,(L+2));
    pg = Data(4,(L+2));
    Declanation = Declanation_list(L);
 for Tilt = -10:80
  for Azimuth = -60:60
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
   aa = a - (hd/h);    

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
   R = D + ((hd/h) * (.5*(1+cosd(Tilt)))) + (pg*(.5*(1-cosd(Tilt))));
   Ht = (h * R);
   Htilt((Tilt+11),(Azimuth+61)) = Ht;
   Alpha((Tilt+11),(Azimuth+61)) = Azimuth;
   Beta((Tilt+11),(Azimuth+61)) = Tilt;
   Tc((Tilt+11),(Azimuth+61)) = Ta + [(NOCT - 20)/800] * Ht;
   Pdc((Tilt+11),(Azimuth+61)) = 0.16 * Ht * (1+((-.4/100)*(Tc((Tilt+11),(Azimuth+61))-25)));
  end
 end
Htiltf = [Htiltf Htilt];
Htiltf(imag(Htiltf) ~= 0) = 0;
Tcf = [Tcf Tc];
Tcf(imag(Tcf) ~= 0) = 0;
Pdcf = [Pdcf Pdc];
Pdcf(imag(Pdcf) ~= 0) = 0;
end 
Janh = Htiltf(1:91 , 1:121);     Ho(1)  = max(Janh(:)); [r1 , c1] = find(Janh == Ho(1),  1, 'first');  To(1) = r1 - 11;   Ao(1) = c1 - 61;                
Febh = Htiltf(1:91, 182:242);    Ho(2)  = max(Febh(:)); [r2 , c2] = find(Febh == Ho(2),  1, 'first');   To(2) = r2 - 11;   Ao(2) = c2 - 61;
Marh = Htiltf(1:91 , 243:363);   Ho(3)  = max(Marh(:)); [r3 , c3] = find(Marh == Ho(3),  1, 'first');  To(3) = r3 - 11;   Ao(3) = c3 - 61;
Aprh = Htiltf(1:91 , 364:484);   Ho(4)  = max(Aprh(:)); [r4 , c4] = find(Aprh == Ho(4),  1, 'first');  To(4) = r4 - 11;   Ao(4) = c4 - 61;
Mayh = Htiltf(1:91 , 485:605);   Ho(5)  = max(Mayh(:)); [r5 , c5] = find(Mayh == Ho(5),  1, 'first');  To(5) = r5 - 11;   Ao(5) = c5 - 61;
Junh = Htiltf(1:91 , 606:726);   Ho(6)  = max(Junh(:)); [r6 , c6] = find(Junh == Ho(6),  1, 'first');  To(6) = r6 - 11;   Ao(6) = c6 - 61;
Julh = Htiltf(1:91 , 727:847);   Ho(7)  = max(Julh(:)); [r7 , c7] = find(Julh == Ho(7),  1, 'first');  To(7) = r7 - 11;   Ao(7) = c7 - 61;
Augh = Htiltf(1:91 , 848:968);   Ho(8)  = max(Augh(:)); [r8 , c8] = find(Augh == Ho(8),  1, 'first');  To(8) = r8 - 11;   Ao(8) = c8 - 61;
Seph = Htiltf(1:91 , 969:1089);  Ho(9)  = max(Seph(:)); [r9 , c9] = find(Seph == Ho(9),  1, 'first');  To(9) = r9 - 11;   Ao(9) = c9 - 61;
Octh = Htiltf(1:91 , 1090:1210); Ho(10) = max(Octh(:)); [r10, c10] = find(Octh == Ho(10),1, 'first'); To(10) = r10 - 11; Ao(10) = c10 - 61;
Novh = Htiltf(1:91 , 1211:1331); Ho(11) = max(Novh(:)); [r11, c11] = find(Novh == Ho(11),1, 'first'); To(11) = r11 - 11; Ao(11) = c11 - 61;
Dech = Htiltf(1:91 , 1332:1452); Ho(12) = max(Dech(:)); [r12, c12] = find(Dech == Ho(12),1, 'first'); To(12) = r12 - 11; Ao(12) = c12 - 61;
for Aoo = 1:12
    if Ao(Aoo) == -60
        Ao(Aoo) = 0;
    end
end 

Janp = Pdcf(1:91 , 1:121);     Po(1)  = max(Janp(:)); [rp1 , cp1] = find(Janp == Po(1),  1, 'first');  Top(1) = rp1 - 11;   Aop(1) = cp1 - 61;                
Febp = Pdcf(1:91 , 122:242);   Po(2)  = max(Febp(:)); [rp2 , cp2] = find(Febp == Po(2),  1, 'first');   Top(2) = rp2 - 11;  Aop(2) = cp2 - 61;
Marp = Pdcf(1:91 , 243:363);   Po(3)  = max(Marp(:)); [rp3 , cp3] = find(Marp == Po(3),  1, 'first');  Top(3) = rp3 - 11;   Aop(3) = cp3 - 61;
Aprp = Pdcf(1:91 , 364:484);   Po(4)  = max(Aprp(:)); [rp4 , cp4] = find(Aprp == Po(4),  1, 'first');  Top(4) = rp4 - 11;   Aop(4) = cp4 - 61;
Mayp = Pdcf(1:91 , 485:605);   Po(5)  = max(Mayp(:)); [rp5 , cp5] = find(Mayp == Po(5),  1, 'first');  Top(5) = rp5 - 11;   Aop(5) = cp5 - 61;
Junp = Pdcf(1:91 , 606:726);   Po(6)  = max(Junp(:)); [rp6 , cp6] = find(Junp == Po(6),  1, 'first');  Top(6) = rp6 - 11;   Aop(6) = cp6 - 61;
Julp = Pdcf(1:91 , 727:847);   Po(7)  = max(Julp(:)); [rp7 , cp7] = find(Julp == Po(7),  1, 'first');  Top(7) = rp7 - 11;   Aop(7) = cp7 - 61;
Augp = Pdcf(1:91 , 848:968);   Po(8)  = max(Augp(:)); [rp8 , cp8] = find(Augp == Po(8),  1, 'first');  Top(8) = rp8 - 11;   Aop(8) = cp8 - 61;
Sepp = Pdcf(1:91 , 969:1089);  Po(9)  = max(Sepp(:)); [rp9 , cp9] = find(Sepp == Po(9),  1, 'first');  Top(9) = rp9 - 11;   Aop(9) = cp9 - 61;
Octp = Pdcf(1:91 , 1090:1210); Po(10) = max(Octp(:)); [rp10, cp10] = find(Octp == Po(10),1, 'first'); Top(10) = rp10 - 11;  Aop(10) = cp10 - 61;
Novp = Pdcf(1:91 , 1211:1331); Po(11) = max(Novp(:)); [rp11, cp11] = find(Novp == Po(11),1, 'first'); Top(11) = rp11 - 11;  Aop(11) = cp11 - 61;
Decp = Pdcf(1:91 , 1332:1452); Po(12) = max(Decp(:)); [rp12, cp12] = find(Decp == Po(12),1, 'first'); Top(12) = rp12 - 11;  Aop(12) = cp12 - 61;

for Aoop = 1:12
    if Aop(Aoop) == -60
        Aop(Aoop) = 0;
    end
end 