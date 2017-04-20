#! /usr/bin/octave -qf

function main(args)
  [polynomial, delta, output, width, height, max] = arguments(args);
  printf("Using polynomial:\np(x) = %s\n", polyout(polynomial, 'x'));
  convergences = newton_basins(polynomial, width, height, delta, max);
  write_output(convergences, output);
endfunction

# Newton basins for root convergence.
function convergences = newton_basins(f, width, height, delta, max)
  step = 0.01;
  fder = polyder(f);

  # This, the image dimensions, represents the function domain
  columns = -width:step:width;
  lines = -height:step:height;

  # Each root will be associated to an integer.
  # Method fails will be associated to the integer zero.
  integer = 0;

  # Roots. Helps keeping track of roots and associated integers. Matrix with format (each row)
  # [root_real_part, root_imaginary_part, associated_integer]
  roots = [[Inf, Inf], integer++];

  # Coordinates (of pixels) and integers.
  # Matrix with each row in the format [pixel_x_coord, pixel_y_coord, associated_integer]
  convergences = [];

  printf("Writing output. This is a good time to grab yourself a coffee...\n");
  for pixel_x_coord = columns
    for pixel_y_coord = lines
      [x_r, y_i] = newton(f, fder, (pixel_x_coord + (pixel_y_coord * i)), delta, max);

      # Method failed. Associate pixel to the integer zero
      if abs(x_r + (y_i * i)) == Inf
        convergences = [convergences; [pixel_x_coord, pixel_y_coord, 0]];
      # Method got an answer
      else
        for row_index = 1:rows(roots)
          #lower_limit = abs(roots(row_index, 1:3)(1) + (roots(row_index, 1:3)(2) * i)) - delta;
          #upper_limit = abs(roots(row_index, 1:3)(1) + (roots(row_index, 1:3)(2) * i)) + delta;
          # If the roots are really close, they must be associated to the same integer
          #if abs(x_r + (y_i * i)) <= upper_limit && abs(x_r + (y_i * i)) >= lower_limit
          if abs(x_r - roots(row_index, 1:3)(1)) <= delta && abs(y_i - roots(row_index, 1:3)(2)) <= delta
            convergences = [convergences; [pixel_x_coord, pixel_y_coord, roots(row_index, 1:3)(3)]];
            break;
          elseif row_index == rows(roots)
            convergences = [convergences; [pixel_x_coord, pixel_y_coord, integer]];
            roots = [roots; [x_r, y_i], integer++];
          endif
        endfor
      endif
    endfor
  endfor
endfunction

# Newton-Raphson method
function [x_r, y_i] = newton(f, fder, x0, delta, max)
  x_now = x0;
  x_r = y_i = gap = Inf;

  # Will do 'max' + 1 iterations. If 'it' reaches 'max', the method has failed to converge.
  for it = 0:max
    derivative_on_point = polyval(fder, x_now);

    # Got a non-zero derivative
    if derivative_on_point != 0
      x_next = x_now - (polyval(f, x_now) / derivative_on_point);
      if gap > abs(x_next - x_now)
        gap = abs(x_next - x_now);
        # We got a good aproximation already?
        if gap <= delta
          x_r = real(x_next);
          y_i = imag(x_next);
          break;
        endif
      endif
      x_now = x_next;
    # Derivative equals zero ===> Ggwp.
    else
      x_r = y_i = Inf;
      break;
    endif
  endfor
endfunction

# Write output file to GNUPlot
function write_output(convergences, output)
  fd = fopen(output, "w");

  convergences = sortrows(convergences, [1, 2]);
  for i = 1:rows(convergences)
    fprintf(fd, "%.8f %.8f %.8f\n", real(convergences(i, 1)), real(convergences(i, 2)), convergences(i, 3));
  endfor
  fclose(fd);
endfunction

# Extract values from CLI
function [polynomial, delta, output, width, heigth, max] = arguments(args)
  # Default values
  delta = 6e-8;
  output = strcat(pwd(), "/outputs/output.txt");
  width = heigth = 3;
  max = 25;

  i = 1; p = 0;
  while i <= length(args)
    if strcmp("-p", args{i}) && (i + 1 <= length(args))
      # Not a number right after the parameter
      string = args{++i}; polynomial = [];
      if length(string) > 1 && (string(1) == '-' || string(1) == '+') && all(isstrprop(string(2:length(string)), "digit"))
        string = string;
      elseif !all(isstrprop(args{i + 1}, "digit"))
        printf("Invalid input format. Was expecting a number after '-p' parameter. Terminating execution.");
        exit;
      endif

      # Build polynomial
      polynomial = [polynomial, str2num(string)]; i++;
      while i <= length(args)
        string = args{i};
        if length(string) > 1 && (string(1) == '-' || string(1) == '+') && all(isstrprop(string(2:length(string)), "digit"))
          polynomial = [polynomial, str2num(string)];
          i++;
        elseif all(isstrprop(string, "digit"))
          polynomial = [polynomial, str2num(string)];
          i++;
        else
          break
        endif
      endwhile
      p = 1;
      continue;
    elseif strcmp("-d", args{i}) && (i + 1 <= length(args))
      if !all(isstrprop(args{i + 1}, "digit"))
        printf("Invalid input format. Was expecting a positive integer number after '-d' parameter. Terminating execution.");
        exit;
      endif
      i++; delta = 1.0 / (10 ** str2num(args{i})); i++;
      continue;
    elseif strcmp("-m", args{i}) && (i + 1 <= length(args))
      if !all(isstrprop(args{i + 1}, "digit"))
        printf("Invalid input format. Was expecting a positive integer number after '-m' parameter. Terminating execution.");
        exit;
      endif
      i++; max = str2num(args{i}); i++;
      continue;
    elseif strcmp("-o", args{i}) && (i + 1 <= length(args))
      i++; output = strcat(pwd(), strcat("/outputs/", args{i})); i++;
      continue;
    elseif strcmp("-w", args{i}) && (i + 1 <= length(args))
      if !all(isstrprop(args{i + 1}, "digit"))
        printf("Invalid input format. Was expecting a positive integer number after '-w' parameter. Terminating execution.");
        exit;
      endif
      i++; width = args{i}; i++;
      continue;
    elseif strcmp("-h", args{i}) && (i + 1 <= length(args))
      if !all(isstrprop(args{i + 1}, "digit"))
        printf("Invalid input format. Was expecting a positive integer number after '-h' parameter. Terminating execution.");
        exit;
      endif
      i++; heigth = args{i}; i++;
      continue;
    else
      printf("Invalid input format. Check out if paramaters are correct. Terminating execution.");
      exit;
    endif
  endwhile

  # Mandatory paramater
  if (p == 0)
    printf("Missing the mandatory paramater '-p'. Terminating execution.");
    exit;
  endif
endfunction

main(argv());
