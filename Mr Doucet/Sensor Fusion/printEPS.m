function printEPS( filename )
    if exist(['Figures/',filename], 'file') ~= 2
        print(['Figures/',filename],'-depsc');
    end
end