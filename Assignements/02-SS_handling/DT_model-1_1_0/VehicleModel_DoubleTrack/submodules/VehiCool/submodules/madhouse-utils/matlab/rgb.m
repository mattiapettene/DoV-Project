%------------------------------------------------------------------------------%
%% Authors
%  - Kristian Jonasson.
%  - Sebastiano Taddei.
%
% Notes
% -----
%  - This file is a modified version of rgb.m by Kristian Jonasson, which can be
%    found at the following link: https://www.mathworks.com/matlabcentral/fileexchange/24497-rgb-triple-of-color-name-version-2.
%------------------------------------------------------------------------------%

% RGB triplets of the 148 CSS3 colours
function out = rgb( colour_name )
  % RGB triplet corresponding to a CSS colour name that can be found in the CSS3
  % specification. They are 148 in total and should be supported by all modern
  % browsers. The colour names are case insensitive, and the function returns
  % the RGB triplet in the range [0,1]. The function can also be called with
  % the argument 'chart' to show a chart of all the colours.
  %
  % Arguments
  % ---------
  %  - colour_name -> string, the name of the colour. If 'chart' is given, a
  %                   chart of all the colours is shown.
  %
  % Outputs
  % -------
  %  - out -> 1x3 array, the RGB triplet.
  %
  % Usage
  % -----
  %  - out = rgb( 'red' ) -> out will be [1 0 0].
  %  - rgb( 'chart' )     -> shows a chart of all the colours.
  %

  % Create persistent variables
  persistent rgb_triplets colour_names

  % Fill the persistent variables only once
  if isempty( rgb_triplets ) || isempty( colour_names )

    % Get the colour name and the corresponding RGB triplet
    [rgb_triplets, colour_names] = getcolors();

    % Convert the colour names to lower case
    colour_names = lower( colour_names );

    % Convert the hex colour names to decimal
    rgb_triplets = reshape( hex2dec(rgb_triplets), [], 3 );

    % Divide the RGB triplets by 255 to get the RGB triplets in the range [0,1]
    rgb_triplets = rgb_triplets / 255;

  end

  % Check if the input argument is 'chart'
  if strcmpi( colour_name, 'chart' )

    showcolors();

  % Otherwise, the input argument is a colour name
  else

    % Find the index of the input colour name in the list of colour names
    colour_index = find( strcmpi( colour_name, colour_names ) );

    % Check if the input colour name is not found
    if isempty( colour_index )
      error( ['Unknown color: ', colour_name] );

    % Otherwise, return the RGB triplet
    else
      out = rgb_triplets(colour_index(1), :);
    end

  end

end

% Show a chart of all 148 CSS3 colours
function showcolors()
  % Show a chart of all the 148 CSS3 colours.
  %
  % Usage
  % -----
  %  - showcolors()
  %

  % Get the colour names
  [~, colour_names] = getcolors();

  % Get the unique colour groups
  colour_groups = {'White', 'Gray', 'Red', 'Pink', 'Orange', 'Yellow', ...
                   'Brown', 'Green', 'Blue', 'Purple', 'Grey'};

  % The rest of this function is not documented, and it is not very readable,
  % but it works. Therefore, I have limited myself to simply fixing the
  % indentation.
  J  = [1, 3, 6, 8, 9, 10, 11];
  fl = lower( colour_groups );
  nl = lower( colour_names );

  for i = 1:length( colour_groups )
    n(i) = strmatch( fl{i}, nl, 'exact' );
  end

  clf;

  p   = get(0,'screensize');
  wh  = 0.6*p(3:4);
  xy0 = p(1:2)+0.5*p(3:4) - wh/2;

  set( gcf, 'position', [xy0, wh] );
  axes( 'position', [0, 0, 1, 1], 'visible', 'off' );
  hold on;

  x = 0;
  N = 0;

  for i = 1:(length( J ) - 1)
    N = max( N, n(J(i + 1)) - n(J(i)) + (J(i + 1) - J(i)) * 1.3 );
  end

  h = 1 / N;
  w = 1 / (length( J ) - 1);
  d = w / 30;

  for col = 1:(length( J ) - 1);

    y = 1 - h;

    for i = J(col):(J(col + 1) - 1)

      t = text( x + w / 2, y + h / 10 , [colour_groups{i} ' colors'] );

      set( t, 'fontw', 'bold', 'vert', 'bot', 'horiz', 'cent', 'fontsize', 10 );

      y = y - h;

      for k = n(i):n(i+1)-1

        c      = rgb( colour_names{k} );
        bright = (c(1) + 2 * c(2) + c(3)) / 4;

        if bright < 0.5
            txtcolor = 'w';
        else
            txtcolor = 'k';
        end

        rectangle( 'position', [x + d, y, w - 2 * d, h], 'facecolor', c );

        t = text( x + w / 2, y + h / 2, colour_names{k}, 'color', txtcolor );

        set( t, 'vert', 'mid', 'horiz', 'cent', 'fontsize', 9 );

        y = y - h;

      end

      y = y - 0.3 * h;

    end

    x = x + w;

  end

end

% Get the colour name and the corresponding hex value
function [hex_values, colour_names] = getcolors()
    % Get the colour names and the corresponding hex values.
    %
    % Outputs
    % -------
    %  - hex_values   -> 148x3 array, the hex values of the colours.
    %  - colour_names -> 148x1 cell array, the names of the colours.
    %
    % Usage
    % -----
    %  - [hex_values, colour_names] = getcolors()
    %

    css = { %White colors
            'FF','FF','FF', 'White'
            'FF','FA','FA', 'Snow'
            'F0','FF','F0', 'Honeydew'
            'F5','FF','FA', 'MintCream'
            'F0','FF','FF', 'Azure'
            'F0','F8','FF', 'AliceBlue'
            'F8','F8','FF', 'GhostWhite'
            'F5','F5','F5', 'WhiteSmoke'
            'FF','F5','EE', 'Seashell'
            'F5','F5','DC', 'Beige'
            'FD','F5','E6', 'OldLace'
            'FF','FA','F0', 'FloralWhite'
            'FF','FF','F0', 'Ivory'
            'FA','EB','D7', 'AntiqueWhite'
            'FA','F0','E6', 'Linen'
            'FF','F0','F5', 'LavenderBlush'
            'FF','E4','E1', 'MistyRose'
            %Grey colors'
            '80','80','80', 'Gray'
            'DC','DC','DC', 'Gainsboro'
            'D3','D3','D3', 'LightGray'
            'C0','C0','C0', 'Silver'
            'A9','A9','A9', 'DarkGray'
            '69','69','69', 'DimGray'
            '77','88','99', 'LightSlateGray'
            '70','80','90', 'SlateGray'
            '2F','4F','4F', 'DarkSlateGray'
            '00','00','00', 'Black'
            %Red colors
            'FF','00','00', 'Red'
            'FF','A0','7A', 'LightSalmon'
            'FA','80','72', 'Salmon'
            'E9','96','7A', 'DarkSalmon'
            'F0','80','80', 'LightCoral'
            'CD','5C','5C', 'IndianRed'
            'DC','14','3C', 'Crimson'
            'B2','22','22', 'FireBrick'
            '8B','00','00', 'DarkRed'
            %Pink colors
            'FF','C0','CB', 'Pink'
            'FF','B6','C1', 'LightPink'
            'FF','69','B4', 'HotPink'
            'FF','14','93', 'DeepPink'
            'DB','70','93', 'PaleVioletRed'
            'C7','15','85', 'MediumVioletRed'
            %Orange colors
            'FF','A5','00', 'Orange'
            'FF','8C','00', 'DarkOrange'
            'FF','7F','50', 'Coral'
            'FF','63','47', 'Tomato'
            'FF','45','00', 'OrangeRed'
            %Yellow colors
            'FF','FF','00', 'Yellow'
            'FF','FF','E0', 'LightYellow'
            'FF','FA','CD', 'LemonChiffon'
            'FA','FA','D2', 'LightGoldenrodYellow'
            'FF','EF','D5', 'PapayaWhip'
            'FF','E4','B5', 'Moccasin'
            'FF','DA','B9', 'PeachPuff'
            'EE','E8','AA', 'PaleGoldenrod'
            'F0','E6','8C', 'Khaki'
            'BD','B7','6B', 'DarkKhaki'
            'FF','D7','00', 'Gold'
            %Brown colors
            'A5','2A','2A', 'Brown'
            'FF','F8','DC', 'Cornsilk'
            'FF','EB','CD', 'BlanchedAlmond'
            'FF','E4','C4', 'Bisque'
            'FF','DE','AD', 'NavajoWhite'
            'F5','DE','B3', 'Wheat'
            'DE','B8','87', 'BurlyWood'
            'D2','B4','8C', 'Tan'
            'BC','8F','8F', 'RosyBrown'
            'F4','A4','60', 'SandyBrown'
            'DA','A5','20', 'Goldenrod'
            'B8','86','0B', 'DarkGoldenrod'
            'CD','85','3F', 'Peru'
            'D2','69','1E', 'Chocolate'
            '8B','45','13', 'SaddleBrown'
            'A0','52','2D', 'Sienna'
            '80','00','00', 'Maroon'
            %Green colors
            '00','80','00', 'Green'
            '98','FB','98', 'PaleGreen'
            '90','EE','90', 'LightGreen'
            '9A','CD','32', 'YellowGreen'
            'AD','FF','2F', 'GreenYellow'
            '7F','FF','00', 'Chartreuse'
            '7C','FC','00', 'LawnGreen'
            '00','FF','00', 'Lime'
            '32','CD','32', 'LimeGreen'
            '00','FA','9A', 'MediumSpringGreen'
            '00','FF','7F', 'SpringGreen'
            '66','CD','AA', 'MediumAquamarine'
            '7F','FF','D4', 'Aquamarine'
            '20','B2','AA', 'LightSeaGreen'
            '3C','B3','71', 'MediumSeaGreen'
            '2E','8B','57', 'SeaGreen'
            '8F','BC','8F', 'DarkSeaGreen'
            '22','8B','22', 'ForestGreen'
            '00','64','00', 'DarkGreen'
            '6B','8E','23', 'OliveDrab'
            '80','80','00', 'Olive'
            '55','6B','2F', 'DarkOliveGreen'
            '00','80','80', 'Teal'
            %Blue colors
            '00','00','FF', 'Blue'
            'AD','D8','E6', 'LightBlue'
            'B0','E0','E6', 'PowderBlue'
            'AF','EE','EE', 'PaleTurquoise'
            '40','E0','D0', 'Turquoise'
            '48','D1','CC', 'MediumTurquoise'
            '00','CE','D1', 'DarkTurquoise'
            'E0','FF','FF', 'LightCyan'
            '00','FF','FF', 'Cyan'
            '00','FF','FF', 'Aqua'
            '00','8B','8B', 'DarkCyan'
            '5F','9E','A0', 'CadetBlue'
            'B0','C4','DE', 'LightSteelBlue'
            '46','82','B4', 'SteelBlue'
            '87','CE','FA', 'LightSkyBlue'
            '87','CE','EB', 'SkyBlue'
            '00','BF','FF', 'DeepSkyBlue'
            '1E','90','FF', 'DodgerBlue'
            '64','95','ED', 'CornflowerBlue'
            '41','69','E1', 'RoyalBlue'
            '00','00','CD', 'MediumBlue'
            '00','00','8B', 'DarkBlue'
            '00','00','80', 'Navy'
            '19','19','70', 'MidnightBlue'
            %Purple colors
            '80','00','80', 'Purple'
            'E6','E6','FA', 'Lavender'
            'D8','BF','D8', 'Thistle'
            'DD','A0','DD', 'Plum'
            'EE','82','EE', 'Violet'
            'DA','70','D6', 'Orchid'
            'FF','00','FF', 'Fuchsia'
            'FF','00','FF', 'Magenta'
            'BA','55','D3', 'MediumOrchid'
            '93','70','DB', 'MediumPurple'
            '99','66','CC', 'Amethyst'
            '8A','2B','E2', 'BlueViolet'
            '94','00','D3', 'DarkViolet'
            '99','32','CC', 'DarkOrchid'
            '8B','00','8B', 'DarkMagenta'
            '6A','5A','CD', 'SlateBlue'
            '48','3D','8B', 'DarkSlateBlue'
            '7B','68','EE', 'MediumSlateBlue'
            '4B','00','82', 'Indigo'
            %Gray repeated with spelling grey
            '80','80','80', 'Grey'
            'D3','D3','D3', 'LightGrey'
            'A9','A9','A9', 'DarkGrey'
            '69','69','69', 'DimGrey'
            '77','88','99', 'LightSlateGrey'
            '70','80','90', 'SlateGrey'
            '2F','4F','4F', 'DarkSlateGrey' };

    hex_values   = css(:,1:3);
    colour_names = css(:,4);
end