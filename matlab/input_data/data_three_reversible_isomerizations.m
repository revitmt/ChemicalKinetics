classdef data_three_reversible_isomerizations < cChemicalSystem
%	v1
%   Copyright 2017, Viktor Reshiak, All rights reserved.    
%
%   Purpose
%   =======
%   Description of the chemical system
%   R1-R2)  A <-> B
%   R3-R4)  B <-> C
%   R5-R6)  C <-> D


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
        function obj = data_three_reversible_isomerizations(varargin)
            for i = 1:2:nargin
                option = varargin{i};
                switch option
                    case 'kinetic_law'
                        kin_law = varargin{i+1};
                end
            end
            obj@cChemicalSystem('name','apnum2015_ex3','kinetic_law',kin_law);
            
            obj.add_reaction(' A <-> B ');
            obj.add_reaction(' B <-> C ');
            obj.add_reaction(' C <-> D ');
            
            obj.rates = [ 1.0d5; 1.0d5; 1.0d2; 1.0d2; 1.0d4; 1.0d4 ];
            
            obj.X0 = [1000; 10; 10; 50];
        end
    end
    
end

