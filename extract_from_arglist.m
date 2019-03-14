function [vout value] = extract_from_arglist(vin,key)

vout = vin;
value = [];

for i=1:length(vin)-1
  if ischar(vin{i}) && strcmpi(vin{i},key)
    value = vin{i+1};
    vout(i:i+1) = [];
    return
  end
end

