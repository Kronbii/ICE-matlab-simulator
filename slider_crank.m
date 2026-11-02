function [s, phi3, w3, alpha3, v, acc] = slider_crank(th, r, l, e, w2, alpha2)
% Kinematics of an offset slider-crank at crank angles th (degrees).
% The crank pin is at (r cos th, r sin th) and the piston pin slides on y = -e.
% Returns the piston pin position s (m), rod angle phi3 (degrees), rod angular
% velocity w3 and acceleration alpha3, and piston velocity v and acceleration acc.
  sy = -(e + r*sind(th));                     % l sin(phi3)
  s = r*cosd(th) + sqrt(l^2 - sy.^2);
  phi3 = atan2(sy, s - r*cosd(th))*180/pi;
  c3 = cosd(phi3); s3 = sind(phi3);
  % y loop: r sin(th) + l sin(phi3) = -e, differentiated twice
  w3 = -r*w2*cosd(th)./(l*c3);
  alpha3 = (r*w2^2*sind(th) - r*alpha2*cosd(th) + l*w3.^2.*s3)./(l*c3);
  % x loop: r cos(th) + l cos(phi3) = s, differentiated twice
  v = -r*w2*sind(th) - l*w3.*s3;
  acc = -r*alpha2*sind(th) - r*w2^2*cosd(th) - l*alpha3.*s3 - l*w3.^2.*c3;
end
