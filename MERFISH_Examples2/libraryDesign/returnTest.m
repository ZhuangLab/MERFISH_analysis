function returnTest(obj)

    try
        display('Return function');
        k = get(obj, 'cutPrimersMaxHomologySelf');
        set(obj, 'cutPrimersMaxHomologySelf', k);
        return;
    catch
       display('Catch function'); 
    end
end