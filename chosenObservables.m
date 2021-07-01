% returns chosen variables in specific order given by local variable chosen_order
function chosen_v = chosenObservables(chosen_v)
chosen_order = ["x", "y", "h"];
isChosen = ismember(chosen_order, chosen_v);
chosen_v = chosen_order(isChosen);