classdef Grid_class
   %   GRID_CLASS is class for grid object, which contains information 
   %   about the graph of the grid as well as some methods on the grid
   
   properties(GetAccess = 'public', SetAccess = 'public')
      N; % number of nodes
      M; % number of lines
      G; % number of generators
      rnc; % non-parsed case
      x; % array of conductivities
      f; % array of power flows
      E; % connectivity matrix ( look up real name )
      A; % adjacency matrix
      error='';
      err=0;
   end
   
   methods
      
      function obj = Grid_class(filename) % Class Constructor
         [rc, info] = load(filename); % load rnc - MATPOWER case
         if (info~=0), fprintf('Error loading case "%s"\n',filename); 
                       obj.error='Loading error'; obj.err=1;
         else
            fprintf('Case "%s" was successfuly loaded\n',filename);
            obj.rnc = rc;
            % Parsing the structure we just loaded
            obj.E = rc.branch(:,1:2); % connectivity matrix
            obj.N = Sz.r(rc.bus(:,1)); % number of buses
            obj.G = Sz.r(rc.gen(:,1)); % number of generators
            obj.M = Sz.r(rc.branch(:,1)); % number of lines
         end
      end
      
      function remove_parallel(obj) % Removes parallel lines from the grid
         
      end
      
      function remove_leaves(obj) % Reduces grid to the state without leaves
      
      end
      
      function gather_gens(obj) % if there are gens at the same node - gathers them
      end
      
      function obj = dcopf(obj) % runs Optimal Power Flow from MATPOWER package
         mpopt = mpoption('OUT_ALL',0,'VERBOSE',1);
         fname = '/Users/pekap/Desktop/matlab/OOP_grid/logs/log';
         solvedcase = '/Users/pekap/Desktop/matlab/OOP_grid/logs/opf_case';
         [obj.rnc, success] = runopf(obj.rnc,mpopt,fname,solvedcase);
         [obj.rnc, success] = rundcopf(obj.rnc,mpopt,fname,solvedcase);
         
         if (success==1), fprintf('OPF was successfuly solved\n');
         else fprintf('Error calculating OPF\n');
              obj.error='OPF error'; obj.err=1;
         end
      end
      
      function N_1_analysis(obj)
      
      end
      
      function N_2_analysis(obj)
      
      end
      
   end

end

