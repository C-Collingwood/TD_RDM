function [d] = flipstruct(d)
f = fieldnames(d);
for k = 1:length(f)
af     = f{k};
d.(af) = d.(af)';
end
end