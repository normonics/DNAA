classdef Neuron < handle
    
    properties
        id = 'neuron';
        
        t = 1;
        u = 0; % activation state of a Neuron
        v = 0; % change detection state
        a = 0; % adaptation state
        aFr = 0; % firing rate adaptation state
        
        fu = 0; % interaction state of a Neuron
        fu2 = 0;
        
        h = -10; % resting level of a Neuron
        hV = 0; % resting level of change detection state
        hA = 0;  % resting level of adaptation state
        
        tau = 30; % time constant of a Neuron
        tauV = 60; % time constant change detection
        tauA = 100; % time constant adaptation
        tauAFr = 100; % time constant adaption firing rate
        
        vCoefficient = 0;
        aCoefficient = 0;
        aFrCoefficient = 0;
        
        noiseCoefficient = 7;%7 for v7 and v9; %10 is the sweet spot for v4
        
        interaction = struct('slope', 1, 'midpoint', 0); % parameters for interaction function
        
        interaction2 = struct('slope', 1, 'midpoint', -10); % parameters for interaction function
        
        
        
        input = {}; % input
        
        changeDetectionBypassInput = {}; % use this for e.g. feedback to preactivate change detectors that itseld does not go through change detection
        
    end
    
    methods
        
        function obj = Neuron(dim1, dim2, dim3, dim4)
            
            if nargin == 1
                for i = 1:dim1
                    
                    obj(i) = Neuron;
                    obj(i).id = [num2str(i)];
                    
                end
            end
            
            if nargin == 2
                for i = 1:dim1
                    for  j = 1:dim2
                        
                        obj(i,j) =  Neuron;
                        obj(i,j).id = [num2str(i) num2str(j)];
                        
                    end
                end
            end
            
            if nargin == 3
                for i = 1:dim1
                    for j = 1:dim2
                        for k = 1:dim3
                            
                            obj(i,j,k) = Neuron;
                            obj(i,j,k).id = [num2str(i) num2str(j) num2str(k)];
                            
                        end
                    end
                end
            end
            
            if nargin == 4
                for i = 1:dim1
                    for j = 1:dim2
                        for k = 1:dim3
                            for l = 1:dim4
                                
                                obj(i,j,k,l) = Neuron;
                                obj(i,j,k,l).id = [num2str(i) num2str(j) num2str(k) num2str(l)];
                                
                            end
                        end
                    end
                end
            end
        end % Constructor Method
        
        function [] = init(obj,T) % initiate t to 1 and all state variables to their resting levels

            
            for i = 1:numel(obj)
                
                if nargin == 2
                   obj(i).u = zeros(1,T+1) ;
                   obj(i).v = zeros(1,T+1) ;
                   obj(i).a = zeros(1,T+1) ;
                   obj(i).aFr = zeros(1,T+1) ;
                   obj(i).fu = zeros(1,T+1);
                   obj(i).fu2 = zeros(1,T+1);
                end
                
                
                obj(i).t = 1;
                obj(i).u(1) = obj(i).h;
                obj(i).v(1) = obj(i).hV;
                obj(i).a(1) = obj(i).hA;
                obj(i).aFr(1) = 0;
                
                obj(i).fu(1) = 1/(1 + exp(-obj(i).interaction.slope.*(obj(i).h-obj(i).interaction.midpoint)));
                obj(i).fu2(1) = 1/(1 + exp(-obj(i).interaction2.slope.*(obj(i).h-obj(i).interaction2.midpoint)));
            end
            
            
        end
        
        function outputArray = output(obj,t)
            
            outputArray = zeros(size(obj));
            
            for i = 1:numel(obj)
                
                if nargin == 1 % if no time argument use current time of Neuron
                    t = obj(i).t;
                end
                
                
                outputArray(i) = obj(i).fu(t);
                
            end
            
            
        end
        
        function outputSum = sum(obj,t)
            
            if nargin == 1 % if no time argument use current time of the first Neuron in the array
                t = obj(1).t;
            end
            
            outputSum = sum(output(obj, t));
            
            for i = 1:ndims(output(obj, t))
                
                outputSum = sum(outputSum);
                
            end
            
            % outputSum = sum(output(obj));
            
        end
        
        function outputPlus = plus(a,b)
            
            outputPlus = output(a) + output(b);
            
        end
        
        function [] = addInput(obj,newInput)
            
            for i = 1:numel(obj)
            
                obj(i).input{size(obj(i).input,1)+1,1} = newInput;  
        
            end
            
        end
        
        function [] = addBypassInput(obj,newInput)
            
            for i = 1:numel(obj)
            
                obj(i).changeDetectionBypassInput{size(obj(i).changeDetectionBypassInput,1)+1,1} = newInput;  
        
            end
            
        end
            
        function [] = step(obj) % function that steps a neuron forward in time
            
            for i = 1:numel(obj)
                
                if iscell(obj(i).input) %check if input is cell array
                    
                    tempInput = 0; % reset the temporary input value to 0
                    
                    % calculate temporary input based on input cell array
                    
                    for currentInput = 1:size(obj(i).input,1)
                        
                        tempSubInput = 0; % tempSubInput keeps track of values to multiply for a given multiplicative synapse
                        
                        for currentSubInput = 1:numel(obj(i).input{currentInput})
                            
                            if isfloat(obj(i).input{currentInput}{currentSubInput}) % check if it is a numerical input
                                
                                if numel(obj(i).input{currentInput}{currentSubInput}) > 1 %check if it is a static input or a time-series
                                    
                                    if obj(i).t <= numel(obj(i).input{currentInput}{currentSubInput}) % if it is a time series, either take the value from current t, or 0 if simulation runs longer than stimulus
                                        
                                        tempSubInput(currentSubInput) = obj(i).input{currentInput}{currentSubInput}(obj(i).t);
                                        
                                    else
                                        
                                        tempSubInput(currentSubInput) = 0;
                                        
                                    end
                                else
                                    
                                    tempSubInput(currentSubInput) = obj(i).input{currentInput}{currentSubInput}; % static inputs are delivered as continuous input
                                    
                                end
                            else
                                
                                tempSubInput(currentSubInput) = sum(obj(i).input{currentInput}{currentSubInput}, obj(i).t); % if it is not a float, hopefully it is a neuron, take its output at time t
                                
                            end
                            
                        end
                        
                        
                        tempInput = tempInput + prod(tempSubInput);
                        
                    end
                    
                else
                    tempInput = obj(i).input;
                    
                end
                
                
                
                % calculate inputs that bypass change detection
                
                if iscell(obj(i).changeDetectionBypassInput)
                    
                    tempBypassInput = 0;
                    
                    
                    
                    for currentInput = 1:size(obj(i).changeDetectionBypassInput,1)
                        
                        tempBypassSubInput = 0; % tempSubInput keeps track of values to multiply for a given multiplicative synapse
                        
                        for currentSubInput = 1:numel(obj(i).changeDetectionBypassInput{currentInput})
                            
                            if isfloat(obj(i).changeDetectionBypassInput{currentInput}{currentSubInput}) % check if it is a numerical input
                                
                                if numel(obj(i).changeDetectionBypassInput{currentInput}{currentSubInput}) > 1 %check if it is a static input or a time-series
                                    
                                    if obj(i).t <= numel(obj(i).changeDetectionBypassInput{currentInput}{currentSubInput}) % if it is a time series, either take the value from current t, or 0 if simulation runs longer than stimulus
                                        
                                        tempBypassSubInput(currentSubInput) = obj(i).changeDetectionBypassInput{currentInput}{currentSubInput}(obj(i).t);
                                        
                                    else
                                        
                                        tempBypassSubInput(currentSubInput) = 0;
                                        
                                    end
                                else
                                    
                                    tempBypassSubInput(currentSubInput) = obj(i).changeDetectionBypassInput{currentInput}{currentSubInput}; % static inputs are delivered as continuous input
                                    
                                end
                            else
                                
                                tempBypassSubInput(currentSubInput) = sum(obj(i).changeDetectionBypassInput{currentInput}{currentSubInput}, obj(i).t); % if it is not a float, hopefully it is a neuron, take its output at time t
                                
                            end
                            
                        end
                        
                        
                        tempBypassInput = tempBypassInput + prod(tempBypassSubInput);
                        
                    end
                    
                else
                    tempBypassInput = obj(i).changeDetectionBypassInput;
                    
                end
                
                
                
                
                
                % implement step
                
                obj(i).u(obj(i).t+1) = obj(i).u(obj(i).t)...
                    + 1/obj(i).tau * (-obj(i).u(obj(i).t)...
                    + obj(i).h...
                    - obj(i).vCoefficient * obj(i).v(obj(i).t)...
                    - obj(i).aCoefficient * obj(i).a(obj(i).t)...
                    - obj(i).aFrCoefficient * obj(i).aFr(obj(i).t)...
                    + tempInput...
                    + tempBypassInput...
                    + obj(i).noiseCoefficient*randn(1));
                %
                %
                
                obj(i).v(obj(i).t+1) = obj(i).v(obj(i).t)...
                    + 1/obj(i).tauV * (-obj(i).v(obj(i).t)...
                    + obj(i).hV + tempInput);
                
                
                
                
                obj(i).a(obj(i).t+1) = obj(i).a(obj(i).t)...
                    + 1/obj(i).tauA * (-obj(i).a(obj(i).t)...
                    + obj(i).hA + obj(i).u(obj(i).t));
                
                
                
                obj(i).aFr(obj(i).t+1) = obj(i).aFr(obj(i).t)...
                    + 1/obj(i).tauAFr * (-obj(i).aFr(obj(i).t)...
                    + output(obj(i), obj(i).t));
                
                
                obj(i).fu(obj(i).t+1) = 1/(1 + exp(-obj(i).interaction.slope.*(obj(i).u(obj(i).t+1)-obj(i).interaction.midpoint)));
                
                obj(i).fu2(obj(i).t+1) = 1/(1 + exp(-obj(i).interaction2.slope.*(obj(i).u(obj(i).t+1)-obj(i).interaction2.midpoint)));
                
                obj(i).t = obj(i).t + 1;
                
            end
            
        end
        
        function [] = plot(obj)
            
            
            for i = 1:numel(obj)
                
                hold on
                plot(obj(i).u,'k')
                
                if obj(i).vCoefficient ~= 0
                    
                    %plot(obj(i).vCoefficient*obj(i).v,'r')
                    
                end
                
                if obj(i).aCoefficient ~= 0
                    
                    plot(obj(i).aCoefficient*obj(i).a,'b')
                    
                end
                
                if obj(i).aFrCoefficient ~= 0
                    
                    plot(obj(i).aFrCoefficient*obj(i).aFr,'g')
                    
                end
                
            end
        end
        
        function [] = map(obj, threshold, t) % find activations above threshold in 4d array and map them with arrows
            %           cla
         
        
        if ndims(obj) ==2
        n = size(obj,1);    
           
        
                
        motionIndices = find(output(obj, t)>threshold);
        
            if isempty(motionIndices)
                
            else
               
               for i = 1:length(motionIndices)
                        
                        [mTemp(i,1) mTemp(i,2)] = ind2sub([n n], motionIndices(i));
                        
               end 
               
              
               
               
               mTemp(:,2) = mTemp(:,2) - mTemp(:,1);
               m = [mTemp(:,1) ones(size(mTemp,1),1) mTemp(:,2) zeros(size(mTemp,1),1)];
               
               for i = 1:length(motionIndices)
                        
                        hold on
                        quiver(m(i,1), m(i,2), m(i,3), m(i,4),'b', 'LineWidth',2, 'MaxHeadSize', 0.5)
                        axis([.5 n+.5 .5 1.5])
                        
               end
               
            end
     
        
        
        
        
        end
            
            
            
        if ndims(obj) == 4
            
        n = size(obj,1);
            if nargin == 2  
                
                motionIndices = find(output(obj)>threshold);
                
            elseif nargin == 3
                
                motionIndices = find(output(obj, t)>threshold);
                
                if isempty(motionIndices)
                    
                    
                    % plot(0,0)
                    
                else
                    
                    
                    for i = 1:length(motionIndices)
                        
                        [m(i,1) m(i,2) m(i,3) m(i,4)] = ind2sub([n n n n], motionIndices(i));
                        
                    end
                    
                    m(:,3) = m(:,3) - m(:,1);
                    m(:,4) = m(:,4) - m(:,2);
                    
                    m = [m(:,2) m(:,1) m(:,4) m(:,3)];
                    
                    
                    for i = 1:length(motionIndices)
                        
                        hold on
                        quiver(m(i,1), m(i,2), m(i,3), m(i,4),'b', 'LineWidth',2)
                        axis([.5 n+.5 .5 n+.5])
                        %quiver(m(i,1), m(i,2), m(i,3), m(i,4),'b', 'LineWidth',2)
                    end
                    
                end
                
            end
            
        end
        end
        
        function [] = connectStim(obj, stim, synapticStrength)
           
            for i = 1:size(obj,1)
                for j = 1:size(obj,2)
                    
                    addInput(obj(i,j), {synapticStrength, stim(i,j,:)});
                    
                end
                
            end
            
            
        end
        
        function [] = setvCoefficient(obj, vCoefficient)
           
            for i = 1:numel(obj)
                
                obj(i).vCoefficient = vCoefficient;
                
            end
            
            
        end
        
        function [] = setaFrCoefficient(obj, aFrCoefficient)
           
            for i = 1:numel(obj)
                
                obj(i).aFrCoefficient = aFrCoefficient;
                
            end
            
        end
        
        function [] = settauAFr(obj, tauAFr)
            
             for i = 1:numel(obj)
                
                obj(i).tauAFr = tauAFr;
                
             end
        end
        
        function [] = copyU(obj1, obj2)
            
            
            if size(obj1) ~= size(obj2)
                error('not the same size, cant copy u')
            else
                
                for i = 1:numel(obj1)
                    
                    obj2(i).u(obj1(i).t) = obj1(i).u(obj1(i).t);
                    
                    obj2(i).fu(obj1(i).t) = 1/(1 + exp(-obj2(i).interaction.slope.*(obj2(i).u(obj1(i).t) -obj2(i).interaction.midpoint)));
                    
                end
                
            end
            
        end
        
        function [] = setInteractionMidpoint(obj, midpoint)
           
            for i = 1:numel(obj)
               
                obj(i).interaction.midpoint = midpoint;
                
            end
            
        end
        
        function [] = setH(obj, h)
            
            for i = 1:numel(obj)
               
                obj(i).h = h;
                
            end
        
        end
        
    end
    
    
end