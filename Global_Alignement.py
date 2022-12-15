""" TP3 where we code a global alignment program
    based on Needleman-Wunsch algorithm """

def simple_display(seq_top, seq_left, score):
    """ Do a simple display """
    print("Seq_top:  {}\nSeq_left: {}\nScore: {}".format(seq_top[::-1], seq_left[::-1], score))

def nice_display(seq_top, seq_left, score):
    """ Do a nice display """
    # What will be printed
    to_print = ""

    # Print seq_top
    # For each positions in a sequence
    for i, _ in enumerate(seq_top):
        # Add the corresponding letter (reverse order)
        to_print += seq_top[-(i+1)]
    # Print a next line after printing the seq_top
    to_print += "\n"

    # Print middle line
    # For each positions in a sequence
    for i, _ in enumerate(seq_top):
        # If it is a match between the two sequences
        if seq_top[-(i+1)] == seq_left[-(i+1)]:
            # Print a pipe
            to_print += "|"
        # Otherwise
        else:
            # Print a space
            to_print += " "
    # Print a next line after printing the middle line
    to_print += "\n"

    # Print seq_left
    # For each positions in a sequence
    for i, _ in enumerate(seq_top):
        # Add the corresponding letter (reverse order)
        to_print += seq_left[-(i+1)]

    # Add the score at the end
    to_print += "\nScore: {}\n".format(score)
    # Print everything!
    print(to_print)

class DynamicMatrix:
    """ Class to generate an empty matrix """
    def __init__(self, seq_top, seq_left, match, mismatch, indel):
        # Self sequence 1
        self.seq_left = seq_left
        # Self sequence 2
        self.seq_top = seq_top
        # Score system
        self.match = match
        self.mismatch = mismatch
        self.indel = indel
        # The Matrix is a list
        self.matrix = []
        # For each line +1 (the init line, so on LEFT sequence)
        for i in range(len(self.seq_left) + 1):
            # Add a new list
            self.matrix.append([])
            # For each column +1 (the init column)
            for j in range(len(self.seq_top) + 1):
                # Add a empty Cell to the last line,
                # the list we just create
                self.matrix[-1].append(Cell())

    # self representation for print
    def __repr__(self):
        # What will be returned
        return "Scores:\n{}\nPrev_pos:\n{}\n\n".format(self.print_scores(), self.print_prev_pos())

    # self representation for print
    def print_scores(self):
        """ Output the values of the matrix """
        # What will be returned
        ret_scores = ".  . "
        # Print top_seq
        for i in self.seq_top:
            ret_scores += "  {} ".format(i)
        # New line
        ret_scores += "\n"
        # For each line
        for ind, i in enumerate(self.matrix):
            # Print seq_left
            if ind > 0:
                ret_scores += "{} ".format(self.seq_left[ind-1])
            else:
                ret_scores += ". "
            # For each column
            for j in i:
                # If this cell has no value
                if j.score is None:
                    # Add a dot to the return
                    ret_scores += (" . ")
                # If this cell is not empty
                else:
                    # Add its content to the return
                    tmp_val = str(j.score)
                    if len(tmp_val) == 1:
                        ret_scores += " " + tmp_val + " "
                    if len(tmp_val) == 2:
                        ret_scores += tmp_val + " "
                    if len(tmp_val) == 3:
                        ret_scores += tmp_val
                # Always add a space after the value we add
                ret_scores += " "
            # End of this line, go to next line
            ret_scores += "\n"
        # Return the content of the Matrix
        return ret_scores

    # self representation for print
    def print_prev_pos(self):
        """ Output the values of the matrix """
        # What will be returned
        ret_prev_pos = ".   .  "
        # Print top_seq
        for i in self.seq_top:
            ret_prev_pos += "     {} ".format(i)
        # New line
        ret_prev_pos += "\n"
        # For each line
        for ind, i in enumerate(self.matrix):
            # Print seq_left
            if ind > 0:
                ret_prev_pos += "{} ".format(self.seq_left[ind-1])
            else:
                ret_prev_pos += ".   "
            # For each column
            for j in i:
                # If this cell has no value
                if j.prev_pos is None:
                    # Add a dot to the return
                    ret_prev_pos += (".   ")
                # If this cell is not empty
                else:
                    # Add its content to the return
                    tmp_val = str(j.prev_pos)
                    ret_prev_pos += tmp_val
                # Always add a space after the value we add
                ret_prev_pos += " "
            # End of this line, go to next line
            ret_prev_pos += "\n"
        # Return the content of the Matrix
        return ret_prev_pos

    def initialize(self):
        """ Initialize the matrix, i.e. fill the first line and column """
        # First cell is 0
        self.matrix[0][0].score = 0
        # First line
        for i in range(1, len(self.seq_top) + 1):
            # Init the value: its position time indel value
            self.matrix[0][i].score = self.indel*i
            # Init the prev_pos: previous position in the line
            self.matrix[0][i].prev_pos = [0, i-1]
        # First column
        for i in range(1, len(self.seq_left) + 1):
            # Init the value: its position time indel value
            self.matrix[i][0].score = self.indel*i
            # Init the prev_pos: previous position in the column
            self.matrix[i][0].prev_pos = [i-1, 0]

    def initialize2(self): #Modifiier les scores pour l'initialisation de la matrice
        """ Initialize the matrix, i.e. fill the first line and column """
        # First cell is 0
        self.matrix[0][0].score = 0
        # First line
        for i in range(1, len(self.seq_top) + 1):
            # Init the value: its position time indel value
            self.matrix[0][i].score = self.indel*i
            # Init the prev_pos: previous position in the line
            self.matrix[0][i].prev_pos = [0, i-1]
        # First column
        for i in range(1, len(self.seq_left) + 1):
            # Init the value: its position time indel value
            self.matrix[i][0].score = self.indel*i
            # Init the prev_pos: previous position in the column
            self.matrix[i][0].prev_pos = [i-1, 0]
    def initialize3(self): #Modifiier les scores pour l'initialisation de la matrice
        """ Initialize the matrix, i.e. fill the first line and column """
        # First cell is 0
        self.matrix[0][0].score = 0
        # First line
        for i in range(1, len(self.seq_top) + 1):
            # Init the value: its position time indel value
            self.matrix[0][i].score = 0
            # Init the prev_pos: previous position in the line
            self.matrix[0][i].prev_pos = [0, i-1]
        # First column
        for i in range(1, len(self.seq_left) + 1):
            # Init the value: its position time indel value
            self.matrix[i][0].score = 0
            # Init the prev_pos: previous position in the column
            self.matrix[i][0].prev_pos = [i-1, 0]

    def compare(self, ntd_a, ntd_b):
        """ Compare to nucleotides and return:
            the match value if they are identical,
            mismatch otherwise
        """
        # If the two nucleotides are the same
        if ntd_a == ntd_b:
            # the match value is returned
            return self.match
        # else, the mismatch value is returned
        return self.mismatch

    def fill_matrix2(self):
        """ Fill-up the matrix """
        # For each cell of the matrix
        for i in range(1, len(self.seq_left) + 1):
            for j in range(1, len(self.seq_top) + 1):
                # The score from the top is its value plus the indel score
                top_score = self.matrix[i-1][j].score + self.indel
                # Position of top
                top_prev_pos = [i-1, j]

                # The score from the top-left is its value plus:
                # the match value if the nucleotides are the same,
                # the mismatch value otherwise
                top_left_score = self.matrix[i-1][j-1].score + self.compare(self.seq_left[i-1], self.seq_top[j-1])
                # Position of top-left
                top_left_prev_pos = [i-1, j-1]

                # The score from the left is its value plus the indel score
                left_score = self.matrix[i][j-1].score + self.indel
                # Position of top
                left_prev_pos = [i, j-1]

                # This cell is the max of the three values
                # and the position of this value
                current_cell_score = min(top_score, top_left_score, left_score)
                # Update the current cell in the matrix
                self.matrix[i][j].score = current_cell_score
                # Update the previous position
                # Is the diagonal the best score?
                # (start with diagonal because match/mismatch are better than gaps)
                if top_left_score == current_cell_score:
                    # Update current cell prev_pos
                    self.matrix[i][j].prev_pos = top_left_prev_pos
                # Else, is the top the best score?
                elif top_score == current_cell_score:
                    # Update current cell prev_pos
                    self.matrix[i][j].prev_pos = top_prev_pos
                # Else, is the left the best score?
                elif left_score == current_cell_score:
                    # Update current cell prev_pos
                    self.matrix[i][j].prev_pos = left_prev_pos
                # This should never occurs
                else:
                    # Error message, in case of...
                    print("There is a bug at position {}".format([i, j]))
    def fill_matrix3(self):
        """ Fill-up the matrix """
        # For each cell of the matrix
        for i in range(1, len(self.seq_left) + 1):
            for j in range(1, len(self.seq_top) + 1):
                # The score from the top is its value plus the indel score
                top_score = max(0,self.matrix[i-1][j].score + self.indel)
                # Position of top
                top_prev_pos = [i-1, j]

                # The score from the top-left is its value plus:
                # the match value if the nucleotides are the same,
                # the mismatch value otherwise
                top_left_score = max(0,self.matrix[i-1][j-1].score + self.compare(self.seq_left[i-1], self.seq_top[j-1]))
                # Position of top-left
                top_left_prev_pos = [i-1, j-1]

                # The score from the left is its value plus the indel score
                left_score = max(0,self.matrix[i][j-1].score + self.indel)
                # Position of top
                left_prev_pos = [i, j-1]

                # This cell is the max of the three values
                # and the position of this value
                current_cell_score = max(top_score, top_left_score, left_score)
                # Update the current cell in the matrix
                self.matrix[i][j].score = current_cell_score
                # Update the previous position
                # Is the diagonal the best score?
                # (start with diagonal because match/mismatch are better than gaps)
                if top_left_score == current_cell_score:
                    # Update current cell prev_pos
                    self.matrix[i][j].prev_pos = top_left_prev_pos
                # Else, is the top the best score?
                elif top_score == current_cell_score:
                    # Update current cell prev_pos
                    self.matrix[i][j].prev_pos = top_prev_pos
                # Else, is the left the best score?
                elif left_score == current_cell_score:
                    # Update current cell prev_pos
                    self.matrix[i][j].prev_pos = left_prev_pos
                # This should never occurs
                else:
                    # Error message, in case of...
                    print("There is a bug at position {}".format([i, j]))

    def fill_matrix(self):
        """ Fill-up the matrix """
        # For each cell of the matrix
        for i in range(1, len(self.seq_left) + 1):
            for j in range(1, len(self.seq_top) + 1):
                # The score from the top is its value plus the indel score
                top_score = self.matrix[i-1][j].score + self.indel
                # Position of top
                top_prev_pos = [i-1, j]

                # The score from the top-left is its value plus:
                # the match value if the nucleotides are the same,
                # the mismatch value otherwise
                top_left_score = self.matrix[i-1][j-1].score + self.compare(self.seq_left[i-1], self.seq_top[j-1])
                # Position of top-left
                top_left_prev_pos = [i-1, j-1]

                # The score from the left is its value plus the indel score
                left_score = self.matrix[i][j-1].score + self.indel
                # Position of top
                left_prev_pos = [i, j-1]

                # This cell is the max of the three values
                # and the position of this value
                current_cell_score = max(top_score, top_left_score, left_score)
                # Update the current cell in the matrix
                self.matrix[i][j].score = current_cell_score
                # Update the previous position
                # Is the diagonal the best score?
                # (start with diagonal because match/mismatch are better than gaps)
                if top_left_score == current_cell_score:
                    # Update current cell prev_pos
                    self.matrix[i][j].prev_pos = top_left_prev_pos
                # Else, is the top the best score?
                elif top_score == current_cell_score:
                    # Update current cell prev_pos
                    self.matrix[i][j].prev_pos = top_prev_pos
                # Else, is the left the best score?
                elif left_score == current_cell_score:
                    # Update current cell prev_pos
                    self.matrix[i][j].prev_pos = left_prev_pos
                # This should never occurs
                else:
                    # Error message, in case of...
                    print("There is a bug at position {}".format([i, j]))

    def global_alignment2(self):
        """ Make a global alignment of two sequences """
        # Creation of two list which are going to contain
        # the two sequences aligned reversed
        al_seq_top = ""
        al_seq_left = ""
        # The end of the matrix is where the algorithm is starting
        maxj,maxi=self.max_mat()
        current_pos = [maxj,maxi]
        # Get the score
        score = self.matrix[current_pos[0]][current_pos[1]].score

        # Backtracking in the matrix
        # The previous position
        prev_pos = []
        # While the beginning of the matrix is not reach
        while self.matrix[current_pos[0]][current_pos[1]].score!=0:
            # Get the previous position
            prev_pos = self.matrix[current_pos[0]][current_pos[1]].prev_pos
            # If the previous x position is the same
            # as the current x position it is a insertion in seq_left
            # so we put a gap
            if current_pos[0] == prev_pos[0]:
                al_seq_left += "-"
            # Else, we write the previous character in the sequence
            else:
                al_seq_left += self.seq_left[prev_pos[0]]
            # If the previous y position is the same
            # as the current y position it is a insertion in seq_top
            # so we put a gap
            if current_pos[1] == prev_pos[1]:
                al_seq_top += "-"
            # Else, we write the previous character in the sequence
            else:
                al_seq_top += self.seq_top[prev_pos[1]]
            # Current position is now the previous position
            current_pos = prev_pos
        return al_seq_top, al_seq_left, score
    def global_alignment(self):
        """ Make a global alignment of two sequences """
        # Creation of two list which are going to contain
        # the two sequences aligned reversed
        al_seq_top = ""
        al_seq_left = ""
        # The end of the matrix is where the algorithm is starting
        current_pos = [len(self.seq_left),len(self.seq_top)]
        # Get the score
        score = self.matrix[current_pos[0]][current_pos[1]].score

        # Backtracking in the matrix
        # The previous position
        prev_pos = []
        # While the beginning of the matrix is not reach
        while prev_pos != [0, 0]:
            # Get the previous position
            prev_pos = self.matrix[current_pos[0]][current_pos[1]].prev_pos
            # If the previous x position is the same
            # as the current x position it is a insertion in seq_left
            # so we put a gap
            if current_pos[0] == prev_pos[0]:
                al_seq_left += "-"
            # Else, we write the previous character in the sequence
            else:
                al_seq_left += self.seq_left[prev_pos[0]]
            # If the previous y position is the same
            # as the current y position it is a insertion in seq_top
            # so we put a gap
            if current_pos[1] == prev_pos[1]:
                al_seq_top += "-"
            # Else, we write the previous character in the sequence
            else:
                al_seq_top += self.seq_top[prev_pos[1]]
            # Current position is now the previous position
            current_pos = prev_pos
        return al_seq_top, al_seq_left, score
    
    def max_mat(self):
        max=0 
        maxi=0
        maxj=0
        for j in range(len(self.seq_left)):
            for i in range(len(self.seq_top)):
                score=self.matrix[j][i].score
                if(score>max):
                    max=score
                    maxi=i
                    maxj=j
        return (maxj,maxi)

class Cell:
    """ Class to generate a cell of the matrix """
    def __init__(self, score=None, prev_pos=None):
        # Value of the cell
        self.score = score
        # Position (x, y) of the previous cell
        self.prev_pos = prev_pos
        # Put an empty array for prev_pos default value
        if self.prev_pos is None:
            self.prev_pos = []

    # self representation for print
    def __repr__(self):
        # Return the content of the cell
        return "Val={}\nPrev_Pos={}\n".format(self.score, self.prev_pos)

def main():
    """ The main of TP3"""
    #Les fonction d'indice 2 servent pour la distance de levenshtein
    #Les fonction d'indice 3 servent pour l'algo de smith-waterman
    mat = DynamicMatrix("ACGGCTAT", "ACTGTAG", 0, 1, 1)
    mat.initialize2()
    #mat.matrix[5][2].score = 42
    #mat.matrix[5][2].prev_pos = [2, 4]
    print(mat)
    #mat.matrix[2][4].score = (5)
    mat.fill_matrix2()
    print(mat)
    al_seq_top, al_seq_left, score = mat.global_alignment()
    nice_display(al_seq_top, al_seq_left, score)
    
    
    mat = DynamicMatrix("TTTACGGCTATTCCC", "ACTGTAG", 2, -1, -2)
    mat.initialize3()
    #mat.matrix[5][2].score = 42
    #mat.matrix[5][2].prev_pos = [2, 4]
    print(mat)
    #mat.matrix[2][4].score = (5)
    mat.fill_matrix3()
    print(mat)
    al_seq_top, al_seq_left, score = mat.global_alignment2()
    nice_display(al_seq_top, al_seq_left, score)

# Launch the main
main()
# Exit without error
exit(0)
# Always put one extra return line
