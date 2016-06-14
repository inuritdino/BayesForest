function obj = set_type(obj,varargin)

    if nargin > 1 && ischar(varargin{1})

        obj.type = varargin{1};
    else
        obj.type = input('Type of the tree is: ','s');
    end
end