classdef data_apnum2015_ex3 < cChemicalSystem
%	v1
%   Copyright 2017, Viktor Reshiak, All rights reserved.    
%
%   Purpose
%   =======
%   Description of the chemical system
%   R1-R2)  A + B <-> C
%   R3-R4)  A + C <-> B
%   R5-R6)  B + C <-> A


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties ( GetAccess='public', SetAccess='protected', Hidden=false )
        X0;     % initial condition
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                              object constructor                                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = data_apnum2015_ex3(varargin)
            for i = 1:2:nargin
                option = varargin{i};
                switch option
                    case 'kinetic_law'
                        kin_law = varargin{i+1};
                end
            end
            obj@cChemicalSystem('name','apnum2015_ex3','kinetic_law',kin_law);
            
            obj.add_reaction(' A + B <-> C ');
            obj.add_reaction(' A + C <-> B ');
            obj.add_reaction(' B + C <-> A ');
            
            obj.rates = [ 1.0d3; 1.0d3; 1.0d-5; 1.0d1; 1.0d0; 1.0d6 ];
            
            obj.X0 = [1e3; 1e3; 1e6];
        end
    end
    
end

