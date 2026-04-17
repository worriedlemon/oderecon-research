function str = method2str(obj)
    if isobject(obj)
        str = class(obj);
    else
        str = func2str(obj);
    end
end