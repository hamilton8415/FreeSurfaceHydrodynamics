classdef LinearIncidentWave < handle
    properties (Access = private)
        ObjectHandle % Holds the uint64 representation of the C++ pointer
    end
    
    methods
        % Constructor
        function obj = LinearIncidentWave()
            % Call MEX to instantiate C++ object
            obj.ObjectHandle = LinearIncidentWave_mex('New');
        end

        %SetSeed method
        function SetSeed(obj, val)
            LinearIncidentWave_mex('SetSeed', obj.ObjectHandle, val);
        end

        %SetToMonoChromatic method
        function SetToMonoChromatic(obj, A, T, phase, beta)
            LinearIncidentWave_mex('SetToMonoChromatic', obj.ObjectHandle, A, T, phase, beta);
        end


        % eta method
        function varargout = eta(obj,varargin)
          if(nargin == 4)
             if(nargout == 1)
               varargout{1} = LinearIncidentWave_mex('eta', obj.ObjectHandle, varargin{1}, varargin{2}, varargin{3});
             end
             if(nargout == 3)
               [varargout{1} varargout{2} varargout{3}] = LinearIncidentWave_mex('eta', obj.ObjectHandle, varargin{1}, varargin{2}, varargin{3});
             end
          end
          if(nargin == 5)
             val = LinearIncidentWave_mex('eta', obj.ObjectHandle, varargin{1}, varargin{2}, varargin{3},varargin{4});
          end
        end

        % MajorVersionNumber method
        function val = MajorVersionNumber(obj)
            val = LinearIncidentWave_mex('MajorVersionNumber', obj.ObjectHandle);
        end

        % MinorVersionNumber method
        function val = MinorVersionNumber(obj)
            val = LinearIncidentWave_mex('MinorVersionNumber', obj.ObjectHandle);
        end

        % PatchVersionNumber method
        function val = PatchVersionNumber(obj)
            val = LinearIncidentWave_mex('PatchVersionNumber', obj.ObjectHandle);
        end 
        
        % Destructor (triggers automatically when MATLAB clears the object)
        function delete(obj)
            if ~isempty(obj.ObjectHandle)
                LinearIncidentWave_mex('Delete', obj.ObjectHandle);
                obj.ObjectHandle = [];
            end
        end
    end
end
