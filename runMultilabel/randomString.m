function string = randomString(stringLength)

symbols = ['a':'z', '0':'9'];
nums = randi(numel(symbols),[1 stringLength]);
string = symbols (nums);

end
