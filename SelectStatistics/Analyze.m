load battery.dsge;
outfile = "selected.dsge";
s = battery;
s = sortbyc(s,1);
s = s(:,2:end);

% see best 5
for i = 1:5
    b = 1:columns(s);
    b = b(:,logical(s(i,:)));
    b
    if i == 1
        save("-ascii", outfile, "b");
    endif        
endfor    

a = mean(s)';
a = [a [1:rows(a)]'];
a = sortbyc(a,1);
a
