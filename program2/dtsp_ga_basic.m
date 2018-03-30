function [ globalMin, optEngy, optRouteT, optRouteD] = dtsp_ga_basic(nStops, popSize, numIter, xy, range, speed  )

% If null arguments
  showprogress=true;
if nargin < 6   % if input variables less than 6, use defualts
  nStops=50;  popSize=400; numIter=500;   speed=10; range=15;     
  xy = 25*rand([nStops, 2]);
  %[nStops, ~] = size(xy);
  showprogress=true;
end

% initialize variables to integers/standard size, create distance matrix
    nPoints = size(xy,1); 
    popSize     =  5*floor(popSize/5);
    numIter     =  max(1,round(real(numIter(1))));
    meshg = meshgrid(1:nPoints);
    dmat = reshape(sqrt(sum((xy(meshg,:)-xy(meshg',:)).^2,2)),nPoints,nPoints);
 
% Initialize the Population
    [n, ~]=  size(xy);
    pop = zeros(popSize,n);
    pop(1,:) = (1:n);
    for k = 2:popSize
        pop(k,:) = randperm(n);
    end
        
 % Run the GA
    globalMin = Inf;
    % assume km 
    truck_energy_km = 10.5e6;
    drone_energy_km =  4.6e4;
    totalDist = zeros(1,popSize);
    totalEngy = zeros(1,popSize);
    truck_route = zeros(popSize, n+1);
    drone_route = zeros(popSize, n+1);
    distHistory = zeros(1,numIter);
    tmpPop = zeros(5,n);
    newPop = zeros(popSize,n);

 for iter = 1:numIter
        % Evaluate Each Population Member (Calculate Total Distance)
        for p = 1:popSize
            k=1; d=0; eng=0;
            while k <= n-1   % dist cities 1 to n
                case_ = 1; % truck only
                  % can one drone be assigned? k+2    
                  if k<n-2 
                        dd1 =  dmat( pop(p,k),     pop(p,k+1) ) + ...
                               dmat( pop(p,k+1  ), pop(p,k+2));
                           if dd1 < range
                               case_ = 2; %assign 2 drones
                           end
                  %else assign truck only    
                  end
             switch case_
                     case 1 % truck only carries drone
                     dt =  dmat( pop(p,k), pop(p,k+1));
                     dd1=0;        % drone rides along, drone distance is zero
                     % truck route carries drone
                     truck_route(p,[k k+1])=[pop(p,k) pop(p,k+1) ];
                     % drones ride inside truck
                     drone_route(p, [k k+1])=[pop(p,k) pop(p,k+1) ];
                     k=k+1; % increment route id by 1
                   case 2 % truck-1-drone, carries other drone
                     dt =  dmat( pop(p,k), pop(p,k+2));
                     % drone is lanched, delivers, recovered
                     truck_route(  p,[k k+1 k+2]) = [pop(p,k) pop(p,k) pop(p,k+2)];
                     % done 1 takes separate path in operation parallel
                     drone_route(p,[k k+1 k+2]) = [pop(p,k) pop(p,k+1) pop(p,k+2)];
                     k=k+2;  % increment route id by 2                          
             end
                d = d + max( [ dt dd1/speed ]); % max of two distance convert to time
                eng = eng + dt*truck_energy_km + dd1*drone_energy_km; % calc energy               
            end 
            % wrap route into circuit back to depot
            if k==n-1
              d1 =    dmat( pop(p,n-1),   pop(p,n)) + ...;
                      dmat( pop(p,n),     pop(p,1));
                truck_route( p,[n-1 n])=[pop(p,n-1) pop(n)];
                drone_route( p,[n-1 n])=[pop(p,n-1) pop(n)];      
            else
              d1 =    dmat( pop(p,n),     pop(p,1));    
            end
            totalDist(p) = d   + d1;                     % total time   
            totalEngy(p) = eng + (d1*truck_energy_km);   % total energy
        end
        
        % Find the Best Route in the Population
        [minDist,index] = min(totalDist);
        distHistory(iter) = minDist;
        if minDist < globalMin
            globalMin = minDist;
            optRouteT = truck_route(index,:);
            optRouteD  = drone_route(index,:);
            optEngy  = totalEngy(index);
            optRouteT=optRouteT(optRouteT>0);
            optRouteD=optRouteD(optRouteD>0);
            if showprogress
                nT = length(optRouteT); 
                nD  = length(optRouteD);
                rtet  = optRouteT([1:nT 1] );   % truck route
                rted  = optRouteD([1:nD 1] );   % drone route
                plot(xy(rted,1),  xy(rted,2), 'k--'); hold on;               
                plot(xy(rtet,1),  xy(rtet,2),'ks-');   
                plot(xy(:,1), xy(:,2),'k.'); 
                xlabel('x-coordinate (km)');
                ylabel('y-coordinate (km)');
                legend('drone','truck','stop', ...
                'Location','bestoutside','Orientation','horizontal')
                 title(sprintf('Truck-1-drone time %1.1f energy  %1.1f MJ',minDist, optEngy/1.0e6)); 
%                 energy = %1.1f MJ',minDist,optEngy/1.0e6));
                hold off;   
                drawnow;
            end  
        end
        
        % Genetic Algorithm Operators
        randomOrder = randperm(popSize);

        for p = 5:5:popSize
                % basically a random sampling in matrix format with a 
            rtes = pop(randomOrder(p-4:p),:); 
            dists = totalDist(randomOrder(p-4:p));
                % what are the min distances?
            [~,idx] = min(dists); 
                % what is the best route
            bestOf5Route = rtes(idx,:);
                % randomly select two route insertion points and sort
            routeInsertionPoints = sort(ceil(n*rand(1,2)));
                I = routeInsertionPoints(1);
                J = routeInsertionPoints(2);
            for k = 1:5 % Mutate the Best row (dist) to get Three New Routes and orig.
                % a small matrix of 4 rows of best time
                tmpPop(k,:) = bestOf5Route;
                switch k
                       % flip two of the cities and cities between
                    case 2 % Flip
                        tmpPop(k,I:J) = tmpPop(k,J:-1:I);
                    case 3 % Swap
                        tmpPop(k,[I J]) = tmpPop(k,[J I]);
                    case 4 % Slide segment 
                        tmpPop(k,I:J) = tmpPop(k,[I+1:J I]);     
                    case 5 % increment sequence one space 
                        tmpPop(k,1:end) = tmpPop(k,[2:end 1]);
                    otherwise % Do Nothing
                end
            end
             % using the original population, create a new population
            newPop(p-4:p,:) = tmpPop;
        end
        pop = newPop;      

 end

end
    