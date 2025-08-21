

function counterStr(i)

    if i > 1
        for blahBlah = 0:log10(i-1)
            fprintf('\b');
        end
    end
    fprintf('%s',num2str(i));
    
end