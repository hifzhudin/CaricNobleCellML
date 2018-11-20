
function [VOI, STATES, ALGEBRAIC, CONSTANTS] = mainFunction()
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
   [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel();
end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".  
    algebraicVariableCount =7;
end
% There are a total of 3 entries in each of the rate and state variable arrays.
% There are a total of 20 entries in the constant variable array.
%

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel()
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;

    % Set timespan to solve over 
    tspan = [0, 5000];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-006, 'AbsTol', 1e-006, 'MaxStep', 1);

    % Solve model with ODE solver
    [VOI, STATES] = ode15s(@(VOI, STATES)computeRates(VOI, STATES, CONSTANTS), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI);

    % Plot state variables against variable of integration
    [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends();
    figure();
    plot(VOI, STATES);
    xlabel(LEGEND_VOI);
    l = legend(LEGEND_STATES);
    set(l,'Interpreter','none');
end

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('time in component main (ms)');
    LEGEND_STATES(:,1) = strpad('V in component main (mV)');
    LEGEND_STATES(:,2) = strpad('h in component main (dimensionless)');
    LEGEND_STATES(:,3) = strpad('n in component main (dimensionless)');
    LEGEND_ALGEBRAIC(:,7) = strpad('i_Na in component main (dimensionless)');
    LEGEND_ALGEBRAIC(:,6) = strpad('i_K in component main (dimensionless)');
    LEGEND_CONSTANTS(:,1) = strpad('G_Na in component main (dimensionless)');
    LEGEND_CONSTANTS(:,2) = strpad('E_Na in component main (dimensionless)');
    LEGEND_CONSTANTS(:,3) = strpad('East in component main (dimensionless)');
    LEGEND_CONSTANTS(:,4) = strpad('Edag in component main (dimensionless)');
    LEGEND_CONSTANTS(:,5) = strpad('g_2 in component main (dimensionless)');
    LEGEND_CONSTANTS(:,6) = strpad('F_n in component main (dimensionless)');
    LEGEND_CONSTANTS(:,7) = strpad('k_1 in component main (dimensionless)');
    LEGEND_CONSTANTS(:,8) = strpad('k_2 in component main (dimensionless)');
    LEGEND_CONSTANTS(:,9) = strpad('k_3 in component main (dimensionless)');
    LEGEND_CONSTANTS(:,10) = strpad('E_1 in component main (dimensionless)');
    LEGEND_CONSTANTS(:,11) = strpad('E_2 in component main (dimensionless)');
    LEGEND_CONSTANTS(:,12) = strpad('E_3 in component main (dimensionless)');
    LEGEND_CONSTANTS(:,13) = strpad('F_h in component main (dimensionless)');
    LEGEND_CONSTANTS(:,14) = strpad('epsilon in component main (dimensionless)');
    LEGEND_CONSTANTS(:,15) = strpad('epsilon_2 in component main (dimensionless)');
    LEGEND_ALGEBRAIC(:,4) = strpad('G in component main (dimensionless)');
    LEGEND_ALGEBRAIC(:,1) = strpad('H_1 in component main (dimensionless)');
    LEGEND_ALGEBRAIC(:,2) = strpad('H_2 in component main (dimensionless)');
    LEGEND_ALGEBRAIC(:,3) = strpad('H_3 in component main (dimensionless)');
    LEGEND_ALGEBRAIC(:,5) = strpad('Istim in component main (dimensionless)');
    LEGEND_CONSTANTS(:,16) = strpad('IstimStart in component main (dimensionless)');
    LEGEND_CONSTANTS(:,17) = strpad('IstimEnd in component main (dimensionless)');
    LEGEND_CONSTANTS(:,18) = strpad('IstimAmplitude in component main (dimensionless)');
    LEGEND_CONSTANTS(:,19) = strpad('IstimPeriod in component main (dimensionless)');
    LEGEND_CONSTANTS(:,20) = strpad('IstimPulseDuration in component main (dimensionless)');
    LEGEND_RATES(:,3) = strpad('d/dt n in component main (dimensionless)');
    LEGEND_RATES(:,2) = strpad('d/dt h in component main (dimensionless)');
    LEGEND_RATES(:,1) = strpad('d/dt V in component main (mV)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [STATES, CONSTANTS] = initConsts()
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    STATES(:,1) = -93.3333;
    STATES(:,2) = 1.0;
    STATES(:,3) = 0.0;
    CONSTANTS(:,1) = 33.33333;
    CONSTANTS(:,2) = 40.0;
    CONSTANTS(:,3) = -15.0;
    CONSTANTS(:,4) = -80.0;
    CONSTANTS(:,5) = -9.0;
    CONSTANTS(:,6) = 0.0037037;
    CONSTANTS(:,7) = 0.075;
    CONSTANTS(:,8) = 0.04;
    CONSTANTS(:,9) = 0.1;
    CONSTANTS(:,10) = -93.3333;
    CONSTANTS(:,11) = -55.0;
    CONSTANTS(:,12) = 1.0;
    CONSTANTS(:,13) = 0.5;
    CONSTANTS(:,14) = 1.0;
    CONSTANTS(:,15) = 1.0;
    CONSTANTS(:,16) = 0;
    CONSTANTS(:,17) = 50000;
    CONSTANTS(:,18) = 80.0;
    CONSTANTS(:,19) = 1000;
    CONSTANTS(:,20) = 1.0;
    if (isempty(STATES)), warning('Initial values for states not set');, end
end

function [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS)
    global algebraicVariableCount;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
    end
    ALGEBRAIC(:,2) = piecewise({STATES(:,1)<=CONSTANTS(:,4), 1.00000 }, 0.000000);
    RATES(:,2) = ( CONSTANTS(:,13).*(ALGEBRAIC(:,2) - STATES(:,2)))./CONSTANTS(:,14);
    ALGEBRAIC(:,3) = piecewise({STATES(:,1)>=CONSTANTS(:,4), 1.00000 }, 0.000000);
    RATES(:,3) =  CONSTANTS(:,15).*CONSTANTS(:,6).*(ALGEBRAIC(:,3) - STATES(:,3));
    ALGEBRAIC(:,1) = piecewise({STATES(:,1)>=CONSTANTS(:,3), 1.00000 }, 0.000000);
    ALGEBRAIC(:,7) = ( CONSTANTS(:,1).*STATES(:,2).*(CONSTANTS(:,2) - STATES(:,1)).*ALGEBRAIC(:,1))./CONSTANTS(:,14);
    ALGEBRAIC(:,4) = piecewise({STATES(:,1)<CONSTANTS(:,4),  CONSTANTS(:,7).*(CONSTANTS(:,10) - STATES(:,1)) , STATES(:,1)>=CONSTANTS(:,3),  CONSTANTS(:,9).*(CONSTANTS(:,12) - STATES(:,1)) },  CONSTANTS(:,8).*(STATES(:,1) - CONSTANTS(:,11)));
    ALGEBRAIC(:,6) =  CONSTANTS(:,5).*ALGEBRAIC(:,3).*(STATES(:,3) .^ 4.00000)+ALGEBRAIC(:,4);
    ALGEBRAIC(:,5) = piecewise({VOI>=CONSTANTS(:,16)&VOI<=CONSTANTS(:,17)&(VOI - CONSTANTS(:,16)) -  (floor(((VOI - CONSTANTS(:,16))./CONSTANTS(:,19)))).*CONSTANTS(:,19)<=CONSTANTS(:,20), CONSTANTS(:,18) }, 0.000000);
    RATES(:,1) = ALGEBRAIC(:,5)+ALGEBRAIC(:,7)+ALGEBRAIC(:,6);
   RATES = RATES';
end

% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI)
    ALGEBRAIC(:,2) = piecewise({STATES(:,1)<=CONSTANTS(:,4), 1.00000 }, 0.000000);
    ALGEBRAIC(:,3) = piecewise({STATES(:,1)>=CONSTANTS(:,4), 1.00000 }, 0.000000);
    ALGEBRAIC(:,1) = piecewise({STATES(:,1)>=CONSTANTS(:,3), 1.00000 }, 0.000000);
    ALGEBRAIC(:,7) = ( CONSTANTS(:,1).*STATES(:,2).*(CONSTANTS(:,2) - STATES(:,1)).*ALGEBRAIC(:,1))./CONSTANTS(:,14);
    ALGEBRAIC(:,4) = piecewise({STATES(:,1)<CONSTANTS(:,4),  CONSTANTS(:,7).*(CONSTANTS(:,10) - STATES(:,1)) , STATES(:,1)>=CONSTANTS(:,3),  CONSTANTS(:,9).*(CONSTANTS(:,12) - STATES(:,1)) },  CONSTANTS(:,8).*(STATES(:,1) - CONSTANTS(:,11)));
    ALGEBRAIC(:,6) =  CONSTANTS(:,5).*ALGEBRAIC(:,3).*(STATES(:,3) .^ 4.00000)+ALGEBRAIC(:,4);
    ALGEBRAIC(:,5) = piecewise({VOI>=CONSTANTS(:,16)&VOI<=CONSTANTS(:,17)&(VOI - CONSTANTS(:,16)) -  (floor(((VOI - CONSTANTS(:,16))./CONSTANTS(:,19)))).*CONSTANTS(:,19)<=CONSTANTS(:,20), CONSTANTS(:,18) }, 0.000000);
end

% Compute result of a piecewise function
function x = piecewise(cases, default)
    set = [0];
    for i = 1:2:length(cases)
        if (length(cases{i+1}) == 1)
            x(cases{i} & ~set,:) = cases{i+1};
        else
            x(cases{i} & ~set,:) = cases{i+1}(cases{i} & ~set);
        end
        set = set | cases{i};
        if(set), break, end
    end
    if (length(default) == 1)
        x(~set,:) = default;
    else
        x(~set,:) = default(~set);
    end
end

% Pad out or shorten strings to a set length
function strout = strpad(strin)
    req_length = 160;
    insize = size(strin,2);
    if insize > req_length
        strout = strin(1:req_length);
    else
        strout = [strin, blanks(req_length - insize)];
    end
end


