function erased_compiled_function(func)
  % erase compiled function with name 'func'

    if exist([func '.dll'])
      clear [func '.dll']
      delete [func '.dll']
    elseif exist ([func '.mexw32'])
      clear [func '.dll']
      delete [func '.mexw32']
    end
