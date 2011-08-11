function texname = getTexName(ts,i)
    texname = ['$' deblank(ts.tex(i,:)) '$'];