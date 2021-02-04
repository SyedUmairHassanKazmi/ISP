
H1 = H_Jan + H_Feb + H_Mar; H2 = H_Apr + H_May + H_Jun; H3 = H_Jul + H_Aug + H_Sep; H4 = H_Oct + H_Nov + H_Dec;
M1 = max(H1(:)); M2 = max(H2(:)); M3 = max(H3(:)); M4 = max(H4(:));
M1+M2+M3+M4




H1 = H_Jan + H_Feb + H_Dec; H2 = H_Apr + H_May + H_Mar; H3 = H_Jul + H_Aug + H_Jun; H4 = H_Oct + H_Nov + H_Sep;

H1 = H_Jan + H_Nov + H_Dec; H2 = H_Apr + H_Feb + H_Mar; H3 = H_Jul + H_May + H_Jun; H4 = H_Oct + H_Aug + H_Sep;
H1 = H_Jan + H_Nov + H_Dec + H_Apr + H_Feb + H_Mar; H2 = H_Jul + H_May + H_Jun + H_Oct + H_Aug + H_Sep;

HToo = [HTo(1)/31, HTo(2)/28, HTo(3)/31, HTo(4)/30, HTo(5)/31, HTo(6)/30, HTo(7)/31, HTo(8)/31, HTo(9)/30, HTo(10)/31, HTo(11)/30, HTo(12)/31];
Poo = [Po(1)/31, Po(2)/28, Po(3)/31, Po(4)/30, Po(5)/31, Po(6)/30, Po(7)/31, Po(8)/31, Po(9)/30, Po(10)/31, Po(11)/30, Po(12)/31];
