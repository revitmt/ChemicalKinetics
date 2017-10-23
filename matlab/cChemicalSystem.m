classdef cChemicalSystem < handle
%	v1
%   Copyright 2017, Viktor Reshiak, All rights reserved.    
%
%   Purpose
%   =======
%   Description of general chemical system
%   


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties ( GetAccess='public', SetAccess='protected', Hidden=false )
        name;           % user-defined name of the chemical system
        kinetic_law;    % kinetic law ( 'mass action' or 'combinatorial' )
        species;        % labels of the chemical species
        reactions;      % string array with reactions
        num_reactions;  % number of species
        num_species;    % number of reactions
        nu;             % nu = ( nu_out - nu_in ); stoichiometric matrix
    end
    properties ( GetAccess='public', SetAccess='public', Hidden=false )
        rates;          % reaction rates
    end
    properties ( GetAccess='public', SetAccess='protected', Hidden=true )
        nu_in;          % reactant stoichiometric matrix (transposed to optimize calculations)
        rev_reactions;  % indices of reversible reactions
        irrev_reactions;% indices of irreversible reactions
        propensities;	% function handle to the propensities
        propJacobian;	% function handle to the Jacobian of the propensities
        relaxRates;     % function handle to relaxation rates of reactions
        A;              % function handle to matrix of the linearization of the ODE system
        str_prop;   	% string vector with propensities
        str_propJac;	% string array with Jacobian of the propensities
        str_relaxRates; % string vector with relaxation rates of reactions
        str_A;          % string array with matrix of the linearization of the ODE system
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                              object constructor                                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = cChemicalSystem(varargin)
            % default parameter values
            obj.name = 'Chemical system';
            
            for i = 1:2:nargin
                option = varargin{i};
                switch option
                    case 'kinetic_law'
                        obj.kinetic_law = varargin{i+1};
                    case 'name'
                        obj.name = varargin{i+1};
                end
            end
            
            obj.species         = strings(0);
            obj.reactions       = strings(0);
            obj.num_reactions   = 0;
            obj.num_species     = 0;
            obj.rev_reactions   = zeros(0,2);
            obj.irrev_reactions = zeros(0,1);
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                  add reaction to the system and build stoichiometry                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function add_reaction(obj,reaction)
            % propensities/Jacobian must be rebuilt every time new reaction is added
            obj.propensities = @obj.first_call_propensities;
            obj.propJacobian = @obj.first_call_propJacobian;
            
            obj.num_reactions = obj.num_reactions + 1;
            
            reaction = strip(reaction);
            
            % split into reactants and products
            ind = strfind( reaction, '<->' );
            % if reversible
            if numel(ind)
                reversible = true;
                reactants_products = strsplit( reaction, '<->' );
                obj.reactions(obj.num_reactions,1)   = replace(reaction,'<->','->');
                obj.reactions(obj.num_reactions+1,1) = replace(reaction,'<->','<-');
                obj.rev_reactions(end+1,:) = obj.num_reactions:obj.num_reactions+1;
            % if irreversible
            else
                reversible = false;
                reactants_products = strsplit( reaction, '->' );
                obj.reactions(obj.num_reactions,1) = reaction;
                obj.irrev_reactions(end+1) = obj.num_reactions;
            end
            
            % split into individual reactants and products
            reactants = strip( strsplit( reactants_products{1}, '+' ) );
            products  = strip( strsplit( reactants_products{2}, '+' ) );

            % find all new species
            for i = 1:length(reactants)
                coef_species = strsplit( reactants{i}, ' ' );
                if length( coef_species ) == 2
                    ind = find( obj.species == coef_species{2} );
                    if ~numel(ind)
                        obj.species(end+1) = coef_species{2};
                        obj.num_species = obj.num_species + 1;
                    end
                else
                    ind = find( obj.species == coef_species{1} );
                    if ~numel(ind)
                        obj.species(end+1) = coef_species{1};
                        obj.num_species = obj.num_species + 1;
                    end
                end
            end
            for i = 1:length(products)
                coef_species = strsplit( products{i}, ' ' );
                if length( coef_species ) == 2
                    ind = find( obj.species == coef_species{2} );
                    if ~numel(ind)
                        obj.species(end+1) = coef_species{2};
                        obj.num_species = obj.num_species + 1;
                    end
                else
                    ind = find( obj.species == coef_species{1} );
                    if ~numel(ind)
                        obj.species(end+1) = coef_species{1};
                        obj.num_species = obj.num_species + 1;
                    end
                end
            end
            num_new_species = length(obj.species) - size(obj.nu,1);
            obj.nu(end+1:end+num_new_species,end+1)    = 0;
            obj.nu_in(end+1,end+1:end+num_new_species) = 0;
            
            % stoichiometry of reactants
            for i = 1:length(reactants)
                coef_species = strsplit( reactants{i}, ' ' );
                if length( coef_species ) == 2
                    obj.nu_in( obj.num_reactions, obj.species == coef_species{2} ) = str2double( coef_species{1} );
                else
                    obj.nu_in( obj.num_reactions, obj.species == coef_species{1} ) = 1;
                end
            end
            
            % stoichiometry of products
            for i = 1:length(products)
                coef_species = strsplit( products{i}, ' ' );
                if length( coef_species ) == 2
                    obj.nu( obj.species == coef_species{2}, obj.num_reactions ) = str2double( coef_species{1} );
                else
                    obj.nu( obj.species == coef_species{1}, obj.num_reactions ) = 1;
                end
            end

            % update stoichiometric vector
            if reversible
                obj.num_reactions = obj.num_reactions + 1;
                obj.nu_in(obj.num_reactions,:) = obj.nu(:,obj.num_reactions-1)';
                obj.nu(:,obj.num_reactions)    = obj.nu_in(obj.num_reactions-1,:)';
                
                obj.nu(1:obj.num_species,obj.num_reactions-1) = obj.nu(1:obj.num_species,obj.num_reactions-1) - obj.nu_in(obj.num_reactions-1,:)';
            end
            obj.nu(1:obj.num_species,obj.num_reactions) = obj.nu(1:obj.num_species,obj.num_reactions) - obj.nu_in(obj.num_reactions,:)'; 
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %            first call of propensity/Jacobian functions with symbolic toolbox            %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function a = first_call_propensities(obj,X)
            obj.build_propensities();
            a = obj.propensities(X);
        end
        function a = first_call_propJacobian(obj,X)
            obj.build_propensities();
            a = obj.propJacobian(X);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                       construct symbolic propensities/Jacobian                          %
        %                 and build function handles for numerical evaluation                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function symb_propensities(obj)
            % symbolic vectors of reaction rates, propensities and species
            c = sym('c',[obj.num_reactions 1]);
            a = sym('x',[obj.num_reactions 1]);
            x = sym('x',[obj.num_species 1]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % construct symbolic propensity vector
            if strcmp( obj.kinetic_law, 'combinatorial' )
                for r = 1:obj.num_reactions
                    a(r) = sym(c(r));
                    for s = 1:obj.num_species
                        for j = 1:obj.nu_in(r,s)
                            a(r) = a(r) * ( x(s) - j + 1 );
                        end
                    end
                end
            elseif strcmp( obj.kinetic_law, 'mass_action' )
                for r = 1:obj.num_reactions
                    a(r) = sym(c(r));
                    for s = 1:obj.num_species
                        a(r) = a(r) * x(s)^obj.nu_in(r,s);
                    end
                end
            end
            
            % convert symbolic propensity vector to string
            prop = '[';
            for r = 1:obj.num_reactions
                prop = strcat( prop, char(a(r)), ';' );
            end
            prop(end) = ']';
            
            % convert string to function handle
            for s = 1:obj.num_species
                prop = replace(prop, char(x(s)), sprintf('X(%d)',s) );
            end
            for r = 1:obj.num_reactions
                prop = replace(prop, char(c(r)), sprintf('obj.rates(%d)',r) );
            end
            prop = strcat('@(obj,X)',prop);
            prop = str2func(prop);
            obj.propensities = @(X)prop(obj,X);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % construct symbolic propensity Jacobian
            Ja = simplify(jacobian(a,x));
            
            % convert symbolic propensity Jacobian to string
            Jas = '[';
            for r = 1:obj.num_reactions
                for s = 1:obj.num_species
                    Jas = strcat( Jas, char(Ja(r,s)), ',' );
                end
                Jas(end) = ';';
            end
            Jas(end) = ']';
            
            % convert string to function handle
            for s = 1:obj.num_species
                Jas = replace(Jas, char(x(s)), sprintf('X(%d)',s) );
            end
            for r = 1:obj.num_reactions
                Jas = replace(Jas, char(c(r)), sprintf('obj.rates(%d)',r) );
            end
            Jas = strcat('@(obj,X)',Jas);
            Jas = str2func(Jas);
            obj.propJacobian = @(X)Jas(obj,X);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                           construct propensities/Jacobian                               %
        %                 and build function handles for numerical evaluation                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function build_propensities(obj)
            obj.str_prop       = strings(obj.num_reactions,1);
            obj.str_propJac    = strings(obj.num_reactions,obj.num_species);
            obj.str_relaxRates = strings(obj.num_reactions,1);
            obj.str_A          = strings(obj.num_species,obj.num_species);
            
            % construct propensity vector
            if strcmp( obj.kinetic_law, 'combinatorial' )
                for r = 1:obj.num_reactions
                    obj.str_prop(r) = sprintf('obj.rates(%d)',r);
                    for s = 1:obj.num_species
                        if obj.nu_in(r,s) == 1
                            obj.str_prop(r) = strcat( obj.str_prop(r), sprintf(' * X(%d)',s) );
                        elseif obj.nu_in(r,s) > 1
                            obj.str_prop(r) = strcat( obj.str_prop(r), sprintf(' * X(%d)', s) );
                            for j = 2:obj.nu_in(r,s)
                                obj.str_prop(r) = strcat( obj.str_prop(r), sprintf(' * ( X(%d) - %d )', s, j-1) );
                            end
                        end
                    end
                end
            elseif strcmp( obj.kinetic_law, 'mass_action' )
                for r = 1:obj.num_reactions
                    obj.str_prop(r) = sprintf('obj.rates(%d)',r);
                    for s = 1:obj.num_species
                        if obj.nu_in(r,s) == 1
                            obj.str_prop(r) = strcat( obj.str_prop(r), sprintf(' * X(%d)',s) );
                        elseif obj.nu_in(r,s) > 1
                            obj.str_prop(r) = strcat( obj.str_prop(r), sprintf(' * X(%d)^%d', s, obj.nu_in(r,s)) );
                        end
                    end
                end
            end

            % convert symbolic propensity vector to string
            prop = "[" + obj.str_prop(1);
            for r = 2:obj.num_reactions
                prop = prop + ";" + obj.str_prop(r);
            end
            prop = prop + "]";
            
            % convert string to function handle
            prop = strcat('@(obj,X)',char(prop));
            prop = str2func(prop);
            obj.propensities = @(X)prop(obj,X);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % construct propensity Jacobian
            if strcmp( obj.kinetic_law, 'combinatorial' )
                for r = 1:obj.num_reactions
                    for s = 1:obj.num_species
                        obj.str_propJac(r,s) = sprintf('obj.rates(%d)',r);
                        
                        for s1 = [ 1:s-1 s+1:obj.num_species ]
                            if obj.nu_in(r,s1) == 1
                                obj.str_propJac(r,s) = obj.str_propJac(r,s) + sprintf(" * X(%d)",s1);
                            elseif obj.nu_in(r,s1) > 1
                                obj.str_propJac(r,s) = obj.str_propJac(r,s) + sprintf(" * X(%d)",s1);
                                for j = 2:obj.nu_in(r,s1)
                                   obj.str_propJac(r,s) = obj.str_propJac(r,s) + sprintf(" * ( X(%d) - %d )", s1, j-1);
                                end
                            end
                        end
                        
                        switch obj.nu_in(r,s)
                            case 0
                                obj.str_propJac(r,s) = 0;
                            case 1
                            case 2
                                obj.str_propJac(r,s) = obj.str_propJac(r,s) + sprintf(" * ( 2*X(%d) - 1 )", s);
                            case 3
                                obj.str_propJac(r,s) = obj.str_propJac(r,s) + sprintf(" * ( 3*X(%d)^2 - 6*X(%d) + 2 )", s, s);
                            case 4
                                obj.str_propJac(r,s) = obj.str_propJac(r,s) + sprintf(" * ( 4*X(%d)^3 - 18*X(%d)^2 + 22*X(%d) - 6 )", s, s, s);
                            otherwise
                                temp = sprintf("( 1/X(%d)",s);
                                for j = 2:obj.nu_in(r,s)
                                    temp = temp + sprintf("+1/(X(%d)-%d)",s,j-1);
                                end
                                temp = temp + " )";
                                obj.str_propJac(r,s) = obj.str_propJac(r,s) + "*" + temp;
                        end
                    end
                end
            elseif strcmp( obj.kinetic_law, 'mass_action' )
                for r = 1:obj.num_reactions
                    for s = 1:obj.num_species
                        obj.str_propJac(r,s) = sprintf('obj.rates(%d)',r);
                        
                        for s1 = [ 1:s-1 s+1:obj.num_species ]
                            if obj.nu_in(r,s1) == 1
                                obj.str_propJac(r,s) = strcat( obj.str_propJac(r,s), sprintf(' * X(%d)',s1) );
                            elseif obj.nu_in(r,s1) > 1
                                obj.str_propJac(r,s) = strcat( obj.str_propJac(r,s), sprintf(' * X(%d)^%d', s1, obj.nu_in(r,s1)) );
                            end
                        end
                        
                        switch obj.nu_in(r,s)
                            case 0
                                obj.str_propJac(r,s) = 0;
                            case 1
                            case 2
                                obj.str_propJac(r,s) = strcat( obj.str_propJac(r,s), sprintf(' * 2 * X(%d)', s) );
                            otherwise
                                obj.str_propJac(r,s) = strcat( obj.str_propJac(r,s), sprintf(' * %d * X(%d)', obj.nu_in(r,s), s) );
                        end
                    end
                end
            end
            
            % convert symbolic propensity Jacobian to string
            Jas = "[" + obj.str_propJac(1,1);
            for s = 2:obj.num_species
                Jas = Jas + "," + obj.str_propJac(1,s);
            end
            for r = 2:obj.num_reactions
                Jas = Jas + ";" + obj.str_propJac(r,1);
                for s = 2:obj.num_species
                    Jas = Jas + "," + obj.str_propJac(r,s);
                end
            end
            Jas = Jas + "]";

            % convert string to function handle
            Jas = strcat('@(obj,X)',char(Jas));
            Jas = str2func(Jas);
            obj.propJacobian = @(X)Jas(obj,X);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % construct matrix of the linearization
            for s1 = 1:obj.num_species
                for s2 = 1:obj.num_species
                    flag = 0;
                    for r = 1:obj.num_reactions
                        if obj.str_propJac(r,s2)~="0"
                            switch abs(obj.nu(s1,r))
                                case 0
                                case 1
                                    flag = 1;
                                    if obj.nu(s1,r)>0
                                        obj.str_A(s1,s2) = strcat( obj.str_A(s1,s2), sprintf(' + (%s)',obj.str_propJac(r,s2)) );
                                    else
                                        obj.str_A(s1,s2) = strcat( obj.str_A(s1,s2), sprintf(' - (%s)',obj.str_propJac(r,s2)) );
                                    end
                                otherwise
                                    flag = 1;
                                    if obj.nu(s1,r)>0
                                        obj.str_A(s1,s2) = strcat( obj.str_A(s1,s2), sprintf(' + %d * (%s)',obj.nu(s1,r),obj.str_propJac(r,s2)) );
                                    else
                                        obj.str_A(s1,s2) = strcat( obj.str_A(s1,s2), sprintf(' - %d * (%s)',abs(obj.nu(s1,r)),obj.str_propJac(r,s2)) );
                                    end
                            end
                        end
                    end
                    if ~flag
                        obj.str_A(s1,s2) = "0";
                    end
                end
            end
            
            % convert symbolic matrix to string
            As = "[" + obj.str_A(1,1);
            for s = 2:obj.num_species
                As = As + "," + obj.str_A(1,s);
            end
            for r = 2:obj.num_species
                As = As + ";" + obj.str_A(r,1);
                for s = 2:obj.num_species
                    As = As + "," + obj.str_A(r,s);
                end
            end
            As = As + "]";

            % convert string to function handle
            As = strcat('@(obj,X)',char(As));
            As = str2func(As);
            obj.A = @(X)As(obj,X);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % construct relaxation rates
            for r1 = 1:obj.num_reactions
                for s = 1:obj.num_species
                    if obj.str_propJac(r1,s)~="0"
                        switch abs( obj.nu(s,r1) )
                            case 0
                            case 1
                                if obj.nu(s,r1) > 0
                                    obj.str_relaxRates(r1) = obj.str_relaxRates(r1) + sprintf(' - %s',obj.str_propJac(r1,s));
                                else
                                    obj.str_relaxRates(r1) = obj.str_relaxRates(r1) + sprintf(' + %s',obj.str_propJac(r1,s));
                                end
                            otherwise
                                if obj.nu(s,r1) > 0
                                    obj.str_relaxRates(r1) = obj.str_relaxRates(r1) + sprintf(' - %d*%s',obj.nu(s,r1),obj.str_propJac(r1,s));
                                else
                                    obj.str_relaxRates(r1) = obj.str_relaxRates(r1) + sprintf(' + %d*%s',abs(obj.nu(s,r1)),obj.str_propJac(r1,s));
                                end
                        end
                    end
                end
            end
            for i = 1:size(obj.rev_reactions,1)
                r1 = obj.rev_reactions(i,1);
                r2 = obj.rev_reactions(i,2);
                for s = 1:obj.num_species
                    % second reaction
                    if obj.str_propJac(r2,s)~="0"
                        switch abs( obj.nu(s,r2) )
                            case 0
                            case 1
                                if obj.nu(s,r2) > 0
                                    obj.str_relaxRates(r1) = obj.str_relaxRates(r1) + sprintf(' - %s',obj.str_propJac(r2,s));
                                else
                                    obj.str_relaxRates(r1) = obj.str_relaxRates(r1) + sprintf(' + %s',obj.str_propJac(r2,s));
                                end
                            otherwise
                                if obj.nu(s,r2) > 0
                                    obj.str_relaxRates(r1) = obj.str_relaxRates(r1) + sprintf(' - %d*%s',obj.nu(s,r2),obj.str_propJac(r2,s));
                                else
                                    obj.str_relaxRates(r1) = obj.str_relaxRates(r1) + sprintf(' + %d*%s',abs(obj.nu(s,r2)),obj.str_propJac(r2,s));
                                end
                        end
                    end
                    obj.str_relaxRates(r2) = obj.str_relaxRates(r1);
                end
            end
            
            % convert symbolic relaxation rates to string
            lambd = "[" + obj.str_relaxRates(1);
            for r = 2:obj.num_reactions
                lambd = lambd + ";" + obj.str_relaxRates(r);
            end
            lambd = lambd + "]";
            
            % convert string to function handle
            lambd = strcat('@(obj,X)',char(lambd));
            lambd = str2func(lambd);
            obj.relaxRates = @(X)lambd(obj,X);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                   print summary                                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function summary(obj)
            fprintf('\tI. Chemical system summary:\n');
            fprintf('\t---------------------------\n\n');
            
            obj.propensities(ones(obj.num_species,1));
            
            tab = table;
            tab.Reactions = obj.reactions;
            tab.Rates = obj.rates;
            tab.Propensities = replace(obj.str_prop,'obj.rates','c');
            tab.Relaxation_rates = replace( obj.str_relaxRates,'obj.rates','c');
            disp(tab);
            
            % Stoichiometry matrix
            R = cell(1,obj.num_reactions);
            S = cell(obj.num_species,1);
            for r = 1:obj.num_reactions
                R{r} = sprintf('R%d',r);
            end
            for s = 1:obj.num_species
                S{s} = strcat( sprintf('X%d (',s), char(obj.species(s)), ')' );
            end
            tab = array2table(obj.nu, 'VariableNames', R, 'RowNames', S);
            fprintf('\n\tII. Stoichiometry matrix:\n');
            fprintf('\t---------------------------\n\n');
            disp(tab);
            
            % Jacobian of propensities
            R = cell(obj.num_reactions,1);
            S = cell(obj.num_species,1);
            for r = 1:obj.num_reactions
                R{r} = sprintf('R%d',r);
            end
            for s = 1:obj.num_species
                S{s} = strcat( sprintf('X%d',s) );
            end
            tab = array2table(replace(obj.str_propJac,'obj.rates','c'), 'VariableNames', S, 'RowNames', R);
            fprintf('\n\tIII. Jacobian of propensities:\n');
            fprintf('\t--------------------------------\n\n');
            disp(tab);
        end
    end
end