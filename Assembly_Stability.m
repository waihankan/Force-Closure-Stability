% The test inputs are referred to NorthWestern modern robotics course,
% "Stablity of an Assemble Project."

% Test-1
%{
massList1 = [1 2 25 35;
            2 5 66 42];
            
contactList1 =[1 0  0  0 pi/2 0.1;
              1 2 60 60 pi   0.5;
              2 0 60  0 pi/2 0.5;
              2 0 72  0 pi/2 0.5];

% This inputs are expected to be unstable objects.
k1 = isStable(massList1, contactList1)
%}

% Test-2
massList2 = [1  2 25 35;
             2 10 66 42];
          
contactList2 =[1 0  0  0 pi/2 0.5;
               1 2 60 60 pi   0.5;
               2 0 60  0 pi/2 0.5;
               2 0 72  0 pi/2 0.5];
          
k2 = isStable(massList2, contactList2);

%{
% for the three-body arch in Figure 12.27 from textbook
massList3 =    [1, 5, 110,  70;
                2, 5, 260, 130;
                3, 5, 410,  70];
mu3 = 0.429;
contactList3 = [1, 0,   0,   0,     pi/2, mu3;
                1, 0,  80,   0,     pi/2, mu3;
                1, 2, 160, 160, -2.81993, mu3;
                1, 2, 180, 100, -2.81993, mu3;
                3, 2, 360, 160, -0.32166, mu3;
                3, 2, 340, 100, -0.32166, mu3;
                3, 0, 440,   0,     pi/2, mu3;
                3, 0, 520,   0,     pi/2, mu3];
k = isStable(massList3, contactList3)
%}

% Main Function
function k = isStable(massList, contactList)
numBodies = size(massList, 1);
numContacts = size(contactList, 1);
Fcontact = zeros(3*numBodies, 2*numContacts);
Fext = zeros(3*numBodies, 1);

    for i = 1 : size(contactList, 1)
        body1 = contactList(i, 1);
        body2 = contactList(i, 2);
        
        % 3D contact points.
        contactPoint = [contactList(i, 3), contactList(i, 4), 0];
        
        % normal_angle into body1.
        normal_angle = contactList(i, 5);
        
         mu = contactList(i, 6);
        
        % if the body is ground(0), we don't need to calculate the Fcontact
        % since we only need to calculate the Fcontact of the contacted
        % bodies, i.e calculate others except the body0.
        if body1 ~= 0
            
            cone_angle = [normal_angle - atan2(mu, 1), ...
                      normal_angle + atan2(mu, 1)];            
            % create 3D normal arrays
            normal1 = [cos(cone_angle(1)), sin(cone_angle(1)), 0];
            normal2 = [cos(cone_angle(2)), sin(cone_angle(2)), 0];
            
            % find the moment.
            moment1 = cross(contactPoint, normal1);
            moment2 = cross(contactPoint, normal2);
            
            % create wrenches.
            F1 = [moment1(3); normal1(1); normal1(2)];
            F2 = [moment2(3); normal2(1); normal2(2)];
            
            % collects the results, F1 and F2 in F matrix.
            Fcontact((body1 - 1)*3 + 1 : body1 * 3, ...    
              (i - 1)*2 + 1)                       = F1;
            Fcontact((body1 - 1)*3 + 1 : body1 *3, ...
              (i - 1)*2 + 2)                       = F2;
        end
        
        if body2 ~= 0
            
            cone_angle = [normal_angle - atan2(mu, 1) + pi, ...
                      normal_angle + atan2(mu, 1) + pi];
                  
            normal1 = [cos(cone_angle(1)), sin(cone_angle(1)), 0];
            normal2 = [cos(cone_angle(2)), sin(cone_angle(2)), 0];
            
            moment1 = cross(contactPoint, normal1);
            moment2 = cross(contactPoint, normal2);
            
            F1 = [moment1(3); normal1(1); normal1(2)];
            F2 = [moment2(3); normal2(1); normal2(2)];
            
            Fcontact((body2 - 1)*3 + 1 : body2 * 3, ...    
              (i - 1)*2 + 1)                       = F1;
            Fcontact((body2 - 1)*3 + 1 : body2 *3, ...
              (i - 1)*2 + 2)                       = F2;
        end
    end
  
  % find external forces acting on the bodies, i.e gravity in this case.
  for j = 1 : size(massList, 1)
      
      % center of mass in 3D
      CoM = [massList(j, 3); massList(j, 4); 0];
      mass = massList(j, 2);
      weight = 9.81 .* mass .* [cos(-pi/2); sin(-pi/2); 0];
      momentofWeight = cross(CoM, weight);
      Fext((j- 1)*3 + 1 : j*3, 1) = [momentofWeight(3); weight(1:2, :)];           
  end 
  
  % check whether Fcontact is full rank or not.
  if rank(Fcontact) == numBodies * 3
      f = ones(1, 2*numContacts);
      A = -1 .* eye(2*numContacts);
      b = -1 .* ones(1, 2*numContacts);
      Aeq = Fcontact;
      beq = -1 .* Fext;
      [k, ~, exitflag] = linprog(f,A,b,Aeq,beq);
      
      if exitflag == -2      
            disp("'The object is not stable.'");
      else
            disp("'The object is stable.'");
      end
  else
      disp("'The object is not stable, rank of Fcontact is not equal to the dimension.'");
  end   
end            