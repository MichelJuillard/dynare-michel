function printMakeCheckOctaveErrMsg(modfilename, err)
    printf("\n");
    printf("********************************************\n");
    printf("*** DYNARE-TEST-OCTAVE ERROR ENCOUNTERED ***\n");
    printf("********************************************\n");
    printf("  WHILE RUNNING MODFILE: %s\n", modfilename);
    printf("                    MSG: %s\n", err.message);
    if (isfield(err, 'stack'))
        printf("                IN FILE: %s\n", err.stack(1).file);
        printf("            IN FUNCTION: %s\n", err.stack(1).name);
        printf("     ON LINE and COLUMN: %d and %d\n",err.stack(1).line,err.stack(1).column);
    end
    printf("*************************************\n\n\n");
end
