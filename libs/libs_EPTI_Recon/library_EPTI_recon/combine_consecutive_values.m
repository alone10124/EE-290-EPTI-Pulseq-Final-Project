function combined_array = combine_consecutive_values(arr)
  % Combines consecutive identical values in an array.
  %
  % Args:
  %   arr: The input array.
  %
  % Returns:
  %   combined_array: A cell array with consecutive identical values combined.

  combined_array = [];
  current_value = arr(1);
  current_count = 1;

  for i = 2:length(arr)
    if arr(i) == current_value
      current_count = current_count + 1;
    else
      combined_array = [combined_array,current_value];
      current_value = arr(i);
      current_count = 1;
    end
  end

  combined_array = [combined_array,current_value]; % Add the last group
end
