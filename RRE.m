function RRE( A, b, fid, sas, format )
% RRE performs Reduced Row Echelon form and writes the output to LaTeX code
%
% INPUTS:
%   A - an n x m matrix to reduce
%   b - (optional) the augmented n x p matrix
%     - default is []
%   fid - (optional) the file id to write a LaTeX file to 
%       - default is 1 (standard out)
%       - use 
%           fid = fopen( '/PATH/TO/file.tex', 'w' ); % or 
%           fid = fopen( '/PATH/TO/file.tex', 'a' );
%       to create a file id. It is up to you to close it with
%       	fclose( fid );
%   sas - (optional) [show all steps]
%       - false is default
%       - false combines all zeroing of column about a pivot in 1 step
%       - true shows the zeroing of the colomn about the pivot for each row
%   format - (optional) the format for printing doubles
%       - default '%.3f': 3 decimal places
%
%
% OUTPUTS:
%
%
% EXAMPLE USAGE:
%   RRE( magic( 4 ), [ 1 2 3 ]' );
%   f = fopen( 'test.tex', 'w' ); RRE( magic( 3 ), [], f ); fclose( f );
%
% KNOWN ISSUES:
%   - There is always a newline added to the end of each step. This will
%   result in an extra equation number at the last line of flalign. 
%   Manually remove the last \\ or use flalign* to supress numbers.
%
%
% NOTES:
%   - LaTeX requires amsmath package (for flalign)
%
%
% HISTORY:
%   2014-02-07a: Paul Romanczyk
%   - Initial version
%   2014-02-07b: Paul Romanczyk
%   - Fixed a bug with the wrong sign being shown in the step
%   - Fixed a bug where making the pivot 1 would only show up on the first
%   step.
%   - Switched from left to right justification of matricies (It looks
%   better)
%   - Made combinining steps (what was intended) for sas = false
%   - Added formatted doubles
%   - Cleaned up the multiline LaTeX description source
%
%
% PUBLIC REPOSITORY:
%   https://github.com/pavdpr/matlab-latex.git
%  
%
% COPYRIGHT:
%   The MIT License (MIT)
%
%   Copyright (c) 2014 Paul Romanczyk
%
%   Permission is hereby granted, free of charge, to any person obtaining a
%   copy of this software and associated documentation files (the 
%   "Software"), to deal in the Software without restriction, including 
%   without limitation the rights to use, copy, modify, merge, publish, 
%   distribute, sublicense, and/or sell copies of the Software, and to 
%   permit persons to whom the Software is furnished to do so, subject to 
%   the following conditions:
%
%   The above copyright notice and this permission notice shall be included
%   in all copies or substantial portions of the Software.
%
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
%   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%   MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
%   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
%   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
%   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%

% input setup
if nargin < 5
    format = '%.3f';
end

if nargin < 4
    sas = false;
end

if nargin < 3
    fid = 1;
end

if numel( fid ) ~= 1
    fid = 1;
end

if fid < 0
    fid = 1;
end

if nargin < 2
    b = [];
end


% error checking
[ n, m ] = size( A );
if numel( b ) ~= 0
    [ N, ~ ] = size( b );
    if N ~= n
        error( 'RRE:invalidBsize', ...
            'b must have the same number of rows as A' );
    end
end

i = 1;

% some LaTeX stuff
fprintf( fid, '\\allowdisplaybreaks\n' ); % allow LaTeX on multiple pages
fprintf( fid, '\\begin{flalign}\n' );
fprintf( fid, RRELaTeXmat( A, b ) );
fprintf( fid, '\t&&\\mbox{Original Matrix}\\\\\n' );

while i <= min( n, m ); 
    % make sure there is at least 1 non-zero in column i
    if ~RREhasPivot( A, i )
        i = i + 1;
    else
        % reduce column i
        [ A, b ] = RREreduce( A, i, b, fid, sas, format );
        i = i + 1;
    end
end

% LaTeX stuff
fprintf( fid, '\\end{flalign}\n' );

end

function [ A, b ] = RREreduce( A, i, b, fid, sas, format )
% RREreduce: reduces column i
%
% INPUTS:
%   A - an n x m matrix to reduce
%   i - the column to work on reducing
%   b - (optional) the augmented n x p matrix
%     - default is []
%   fid - (optional) the file id to write a LaTeX file to 
%       - default is 1 (standard out)
%       - use 
%           fid = fopen( '/PATH/TO/file.tex', 'w' ); % or 
%           fid = fopen( '/PATH/TO/file.tex', 'a' );
%       to create a file id. It is up to you to close it with
%       	close( fid );
%   sas - (optional) [show all steps]
%       - false is default
%       - false combines all zeroing of column about a pivot in 1 step
%       - true shows the zeroing of the colomn about the pivot for each row
%   format - (optional) the format for printing doubles
%       - default '%.3f': 3 decimal places
% OUTPUTS:
%   A - reduced form of input A
%   b - reduced form of input b
%
  
    if nargin < 3 
        b = [];
    end
    if nargin < 4
        fid = 1;
    end
    if nargin < 5
        % show all steps
        sas = false;
    end
    if nargin < 6
        format = '%.3f';
    end
    
    [ n, ~ ] = size( A );
    
    % get the pivot
    idx = RREgetPivot( A, i );
    
    if idx ~= i
        % swap row idx and pivot
        [ A, b ] = RREswap( A, b, i, idx );
        fprintf( fid, RRELaTeXmat( A, b ) );
        fprintf( fid, '\t&& R_{%d}\\leftrightarrow R_{%d}\\\\\n', i, idx );
    end
    
    if A( i, i ) == 0.0
        % should never get here due to RREhasPivot in RRE
        error( 'RRE:RREreduce:zeroAtPivot', 'Zero at Pivot' );
    end

    % put a one at the pivot
    [ A, b, s ] = RREone( A, b, i ); 
    if s ~= 1.0
        fprintf( fid, RRELaTeXmat( A, b ) );
        fprintf( fid, [ '\t&& R_{%d}=' format '\\cdot R_{%d}\\\\\n' ], ...
            i, s, i );
    end
    
    % put a zero everywhere else
    if ~sas
        msg = '\t&&\\begin{array}{l}\n';
    end
    for j = 1:n
        if ( ( j ~= i ) && ( A( j, i ) ~= 0 ) )
            [ A, b, s ] = RREzero( A, b, i, j );
            if sas
                fprintf( fid, RRELaTeXmat( A, b ) );
                
                % write what we did
                fprintf( fid, '\t&& R_{%d}=R_{%d}', j, j );
                if s < 0
                    fprintf( fid, [ '+', format ], -s );
                else
                    fprintf( fid, [ '-', format ], s );
                end
                fprintf( fid, '\\cdot R_{%d}\\\\\n', i );
            else               
                msg = strcat( msg, sprintf( '\t\tR_{%d}=R_{%d}', j, j ) );
                if s < 0
                    msg = strcat( msg, sprintf( [ '+', format ], -s ) );
                else
                    msg = strcat( msg, sprintf( [ '-', format ], s ) );
                end
                msg = strcat( msg, '\\cdot R_{' );
                msg = strcat( msg, sprintf( '%d', i ) );
                msg = strcat( msg, '}\\\\\n' );
            end 
        end
    end
    if ~sas
        msg = strcat( msg, '\t\\end{array}\\\\\n' );
        fprintf( fid, RRELaTeXmat( A, b ) );
        fprintf( fid, msg );
    end
end


function tf = RREhasPivot( A, i )
% RREhasPivot: reduces column i
%
% INPUTS:
%   A - an n x m matrix to reduce
%   i - the column to work on reducing
%
% OUTPUTS:
%   tf - true of false if column i has a pivot at or below row i
%
    if sum( abs( A( i:end, i ) ) ) == 0
        tf = false;
    else 
        tf = true;
    end
end

function [ A, b, s ] = RREone( A, b, i )
% RREone: puts a 1 at position A(i,i) using elemetary matrix operations
%
% INPUTS:
%   A - an n x m matrix to reduce
%   b - (optional) the augmented n x p matrix
%     - default is []
%   i - the column to work on reducing
%   
% OUTPUTS:
%   A - reduced form of input A
%   b - reduced form of input b
%   s - the scale such that A(:,i)_output = s .* A(:,i)_input
%
% WARNINGS:
%   An error will be thrown if there is a zero at the pivot location
%
    s = A( i, i );
    if s == 0
        error( 'RRE:RREone:zeroAtPivot', [ 'Zero at Pivot A(', ...
            num2str( i ), ', ', num2str( i ), ')' ] );
    end
    s = 1 / s;
    A( i, : ) = A( i, : ) .* s;
    if numel( b ) > 0
        b( i, : ) = b( i, : ) .* s;
    end
end

function [ A, b, s ] = RREzero( A, b, i, j )
% RREzero: puts a 1 at position A(i,i) using elemetary matrix operations
%
% INPUTS:
%   A - an n x m matrix to reduce
%   b - (optional) the augmented n x p matrix
%     - default is []
%   i - the row that contains the pivot
%   j - the row to reduce around the pivot
%   
% OUTPUTS:
%   A - reduced form of input A
%   b - reduced form of input b
%   s - the scale such that 0 = A(j,i)_output = A(j,i) - s * A(i,i)
%
% WARNINGS:
%   - An error will be thrown if there is a zero at the pivot location
%   - An error will be thrown if you try to reduce a row with itself (i==j)
%

    if A( i, i ) == 0
        error( 'RRE:RREzero:zeroAtPivot', [ 'Zero at Pivot A(', ...
            num2str( i ), ', ', num2str( i ), ')' ] );
    elseif i == j
        error( 'RRE:RREzero:selfReduce', [ 'Trying to reduce row', ...
            num2str( i ), 'with itself' ] );
    end
    s = A( j, i ) / A( i, i );
    A( j, : ) = A( j, : ) - s .* A( i, : );
    if numel( b ) > 0
        b( j, : ) = b( j, : ) - s .* b( i, : );
    end
end

function idx = RREgetPivot( A, i )
% RREgetPivot finds the index to make 1 for column i.
%
% INPUTS:
%   A - an n x m matrix to reduce
%   i - the row that contains the pivot
%   
% OUTPUTS:
%   idx -  the row to reduce
%
% NOTES: 
%   - for numerical stability, we choose the largest absolute value
%

    idx = find( abs( A( i:end, i ) ) == max( abs( A( i:end, i ) ) ), ...
        1, 'First' ) + i - 1;
end

function [ A, b ] = RREswap( A, b, i, j )
% RREswap: swap row i with row j
%
% INPUTS:
%   A - an n x m matrix to reduce
%   b - (optional) the augmented n x p matrix
%     - default is []
%   i - a row to swap
%   j - the other row to swap
%   
% OUTPUTS:
%   A - reduced form of input A
%   b - reduced form of input b
%
    if i == j
        return;
    end
    tmp = A( i, : );
    A( i, : ) = A( j, : );
    A( j, : ) = tmp;
    if numel( b ) > 0
        tmp = b( i, : );
        b( i, : ) = b( j, : );
        b( j, : ) = tmp;
    end
end

function LaTeX = RRELaTeXmat( A, b, format )
% RRELaTeXmat: writes a string containing LaTeX for augmented matrix [A|b]
%
% INPUTS:
%   A - an n x m matrix to reduce
%   b - (optional) the augmented n x p matrix
%     - default is []
%   format - (optional) the format for printing doubles
%       - default '%.3f': 3 decimal places
%   
% OUTPUTS:
%   LaTeX - A string containing LaTeX code for matrix [A|b]
%
    if nargin < 3
        format = '%.3f';
    end
    
    if nargin < 2 
        p = 0;
    else
        p = size( b, 2 );
        A = [ A, b ];
    end
    
    [ n, m ] = size( A );
    
    LaTeX = '\t\\left[\\begin{array}{';
    % setup the forrmating of the array
    for i = 1:(m-p)
        LaTeX = strcat( LaTeX, 'r' );
    end
    if p > 0 
        LaTeX = strcat( LaTeX, '|' );
        for i = 1:p
            LaTeX = strcat( LaTeX, 'r' );
        end
    end
    
    LaTeX = strcat( LaTeX, '}\n' );
    
    % write the matrix
    for i = 1:n
        LaTeX = strcat( LaTeX, '\t\t' );
        for j = 1:m
            LaTeX = strcat( LaTeX, sprintf( format, A( i, j ) ) );
            if j ~= m
                LaTeX = strcat( LaTeX, ' & ' );                
            end
        end
        LaTeX = strcat( LaTeX, ' \\\\\n' );
    end
    LaTeX = strcat( LaTeX, '\t\\end{array}\\right]\n' );
end