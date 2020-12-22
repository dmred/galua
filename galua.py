import sys
import math

class GaloisField():
    """ A finite field, or Galois field.
        Args:
            p (int): A prime number, the base of the field.
            n (int): An exponent representing the degree of a field extension.
                     By default n = 1 (prime field).
            coefs (list): A list of integers representing the coefficients of
                          an irreducible primitive polynomial of degree n over
                          GF(p). Default is the empty list for prime fields.
        Attributes:
            p (int): The prime dimension of the field
            n (int): The degree of the field extension 
            dim (int): The full order of the field, :math:`p^n`.
            w (pthRootOfUnity): The :math:`p^{\\text{th}}` root of unity.
            coefs (list): The coefficients of the irreducible polynomial
            elements (list): A list of all FieldElements in this finite field.
            is_sdb (bool): A boolean which tells us whether the elements'
                           expansion coefficients are in the self-dual
                           basis (True) or the polynomial basis (False). 
                           The default is False.
    """
    def __init__(self, p, n = 1, coefs = []):
        # TODO implement check for prime number
        self.p = p

        # Field extension parameter. 
        # Base prime field if n = 1, otherwise the field is GF(p^n)
        if n > 1:
            self.n = n
        elif n != 1:
            print("Error, invalid exponent for field extension.")
            print("Please enter an exponent n of 2 or greater.\n")
            sys.exit()
        else:
            self.n = 1
            self.coefs = []
        
        # Set separate parameter for the field dimension
        self.dim = int(math.pow(p, n))

        # Initialize the pth root of unity
        self.w = pthRootOfUnity(p)

        # Initialize the coefficients for the irreducible polynomial 
        # to do the field extension with. 
        if len(coefs) > 0:
            # We should have n + 1 coefficients for GF(p^n)
            if len(coefs) == n + 1:
                self.coefs = coefs
            else:
                print("Error, invalid number of coefficients in the irreducible polynomial.")
                print("Field of size " + str(self.p) + "^" + str(self.n) + " should have ", end = "")
                print(str(n + 1) + " coefficients in its irreducible polynomial.")
                sys.exit()


        # Generate the actual field elements
        if self.n == 1:
            # Prime case is easy. No field basis, just the numbers from 0 to p,
            # stored as FieldElements.
            self.elements = []
            for i in range(0, p):
                self.elements.append(FieldElement(self.p, self.n, [i]))
        else:
            # Use the irreducible polynomial to generate the field elements
            # They'll be stored in order as a list of coefficients in the polynomial basis
            # e.g. in dim 4, x^2 + x + 1 is the polynomial, use the basis (1, x) and store
            # the elements as:
            # 0 -> [0, 0], 1 -> [1, 0], x -> [1, 0], x^2 = [1, 1]

            # Hold all the coefficients for each element
            # For simplicity, rather than a list of list, represent each field element as a 
            # string of coefficients, i.e. [0, 1, 1] -> "011"  
            field_list = []

            # The polynomial basis contains n elements
            # The first element is always 0
            self.elements = []
            self.elements.append(FieldElement(self.p, self.n, [0]*self.n))
            field_list.append("0," * (self.n - 1) + "0")

            # The next few elements are initial terms in the poly basis (i.e. x, x^2 ...)
            for i in range(1, self.n):
                next_coefs = [0]*(i) + [1] + [0]*(self.n - i - 1) 
                self.elements.append(FieldElement(self.p, self.n, next_coefs))
                field_list.append(",".join([str(x) for x in next_coefs]))

            # For the n^th power of x, we need to use the irreducible polynomial
            nth_coefs = [((-1) * self.coefs[i]) % self.p for i in range(0, self.n)]
            self.elements.append(FieldElement(self.p, self.n, nth_coefs))
            field_list.append(",".join([str(x) for x in nth_coefs]))

            # For the remaining powers, multiply the previous element by primitive element
            for el in range(self.n + 1, self.dim):
                # Shift all coefficients ahead by 1 power of x and take the sum because
                # we know all the previous elements, and will never get anything 
                # with such a high exponent we don't know it's basis coefficients
                next_coefs = [0] + self.elements[el - 1].exp_coefs
                
                # Get a list of the powers whose coefficients aren't 0
                which_to_sum = [self.elements[i] * co for i, co in enumerate(next_coefs) if co != 0]
                sum = self.elements[0]

                for sum_el in which_to_sum:
                    sum = sum + sum_el

                # TODO Make sure that this element is not already in the list - if it is, then
                # we did not use a true primitive polynomial.
                str_rep = ",".join([str(x) for x in sum.exp_coefs])
                if str_rep not in field_list:
                    self.elements.append(sum)
                    field_list.append(str_rep)
                else:
                    raise ValueError("Repeated field element detected; please make sure your irreducible polynomial is primitive.")
                 
            # This is really dumb, but make sure each element holds a copy of the whole
            # list of the field elements. This makes field multiplication way easier.
            for i in range(len(self.elements)):
                (self.elements[i]).field_list = field_list 
                (self.elements[i]).prim_power = i

        # SDB information
        self.is_sdb = False # Have we indicated an sdb?
        self.sdb = [] # The indices of the elements that make up the sdb
        self.sdb_norms = [] # The trace of the sdb squared - usually 1, but
                            # if the sdb is almost sd, then one is not 1.


    def __getitem__(self, idx):
        """ Access specific elements in the finite field.
            Args:
                idx (int): The index of the element to retrieve. For primes 
                    this is the same as the number itself; for power-of-primes
                    it represents the power of the primitive element.
          
            Returns:
                The element at the specified index in the field. 
              None if idx is out of bounds.
        """
        if idx < self.dim and idx >= (-1 * self.dim):
            return self.elements[idx]
        else:
            print("Error, element out of bounds.")


    def __iter__(self):
        """ Make the finite field iterable. 
            Returns:
                An iterator to the field elements.
        """
        return iter(self.elements)


    def to_sdb(self, sdb_element_indices):
        """ Transform the expansions coefficients to the self-dual basis.
            Currently valid only for fields whose orders are powers of 2.
            Args:
                sdb_element_indices (list): The indices of the FieldElements 
                    (as powers of the primitive element) that represent the
                    self-dual basis. e.g. if the self-dual basis is 
                    :math:`\{ \sigma^3, \sigma^5, \sigma^6 \}`, this list
                    would be [3, 5, 6].
            TODO:
                Implement automatic construction of some self-dual basis.
        """

        if self.n == 1:
            print("Cannot take self-dual basis of a prime field.")
            return

        # Make sure that the provided sdb is valid. In qudit cases, we may
        # also be shuffling the elements, so make sure to get the shuffled copy.
        valid_sdb, valid_element_indices, valid_sdb_norms = self.verify_sdb(sdb_element_indices)

        if not valid_sdb:
            print("Invalid self-dual basis provided.")
            return

        if valid_element_indices != sdb_element_indices:
            print("The order of your self-dual basis elements has changed.")
            print("This is due to the presence of a non-1 normalization coefficient.")
            print("New ordering is " + str(valid_element_indices) + ".")

        # Set the sdb 
        self.is_sdb = True
        self.sdb = valid_element_indices
        self.sdb_norms = valid_sdb_norms

        first_norm_inverse = 1  

        if self.p % 2 == 1 and self.n % 2 == 0: # If p odd, n even
            for i in range(self.p):
                if (self.sdb_norms[0] * i) % self.p == 1:
                    first_norm_inverse = i

        # If all goes well, we can start computing the coefficients
        # in terms of the new elements by using the trace and multiplication
        # functions.
        sdb_els = [self.elements[self.sdb[i]] for i in range(0, self.n)]
        sdb_field_list = []
        for element in self.elements:
            sdb_coefs = [] # Expansion coefficients in the sdb

            for i in range(len(sdb_els)):
                if i == 0:
                    sdb_coefs.append((first_norm_inverse * tr(element * sdb_els[i])) % self.p)
                else:
                    sdb_coefs.append(tr(element * sdb_els[i]))


            sdb_field_list.append(",".join([str(x) for x in sdb_coefs]))

            element.is_sdb = True
            element.sdb_coefs = sdb_coefs

        # Finally, give every element a copy of the sdb coefficients
        for element in self.elements:
            element.sdb_field_list = sdb_field_list
    


    def verify_sdb(self, sdb_element_indices):
        """ Verify if a set of elements form a proper self-dual normal basis.
            For qubit systems, a proper self-dual basis always exists. There are
            two properties to check for this:
              * The trace of each basis element squared is 1.
              * The trace of each basis element multiplied by every other is 
                0 (orthogonality). 
            In qudit cases, the first condition needs to be modified; there
            will always be some basis element where the trace of the square is *not*
            1, but rather some other element in the prime field. To this end, we will
            keep a list of the 'normalization' coefficients for the almost sdb. We will
            also reorder the almost sdb in this case so that the non-1 element is first.
 
            Returns:
                A triple containing the following values:
                - True if above conditions are satisfied, false if not. 
                - The sdb element indices, the ordering of which may change if 
                  the basis is not perfectly self-dual. None if cond 1 is false.
                - The normalizations of the sdb elements. A list of 1s if the
                  basis is perfectly self-dual, or a positive coefficient plus
                  the rest of the list 1s if almost self-dual. None if cond 1
                  is not satisfied.
        """

        if len(sdb_element_indices) != self.n:
            print("Error, incorrect number of elements in proposed basis.")
            return False, None, None

        sdb_norms = []

        if self.p == 2: # Qubit case
            for i in range(0, self.n):
                for j in range(i, self.n): # Don't double compute things
                    trace_result = tr(self.elements[sdb_element_indices[i]] * self.elements[sdb_element_indices[j]])

                    if i == j: # Same element, should have trace 1
                        if trace_result != 1:
                            return False, None, None
                    else: # Different elements should be orthogonal and have trace 0
                        if trace_result != 0:
                            return False, None, None

            # If successful, set the orthogonality coefficients to 1
            sdb_norms = [1] * self.n

        else: # Qudit case
            for i in range(0, self.n):
                for j in range(0, self.n):
                    trace_result = tr(self.elements[sdb_element_indices[i]] * self.elements[sdb_element_indices[j]])
                    
                    if i == j: # Square the element and trace it
                        # Just needs to be in the prime_field
                        if trace_result < 0 or trace_result >= self.p:
                            return False, None, None
                           
                        # Only one element can have a non-1 normalization 
                        if trace_result == 1: 
                            sdb_norms.append(trace_result)
                        else:
                            non1 = [x for x in sdb_norms if x != 1]
                            if len(non1) >= 1:
                                print("Error, more than one element has a non-one normalization coefficient.")
                                print("Self-dual basis is invalid.")
                                return False, None, None
                            else:
                                sdb_norms.append(trace_result)
                    else: # Different elements must be trace-orthogonal
                        if trace_result != 0:
                            return False, None, None

            # For power of primes, the self-dual basis **might** be real (e.g. dim 27).
            # It's possible that all normalizations are 1, so check this, and carry on if true.
            if sdb_norms.count(1) != len(sdb_norms):
                # If the sdb so far has been okay, let's reshuffle it so the element
                # with coefficient > 1 is at the beginning. I'm honestly not sure why we 
                # do this, but this is what Andrei said to do in our LS paper.
                # Get the index of the non-1 element. Thanks to 
                # http://stackoverflow.com/questions/4111412/how-do-i-get-a-list-of-indices-of-non-zero-elements-in-a-list
                non1 = [i for i, e in enumerate(sdb_norms) if e != 1][0] 

                shuffled_sdb = [sdb_element_indices[non1]] + sdb_element_indices[:non1] + \
                                    sdb_element_indices[non1 + 1:]
                shuffled_norms = [sdb_norms[non1]] + sdb_norms[:non1] + sdb_norms[non1 + 1:]

                sdb_element_indices = shuffled_sdb
                sdb_norms = shuffled_norms
                              
        # If we made it this far, we're golden.
        return True, sdb_element_indices, sdb_norms


    def compute_sdb(self):
        """ Compute a self-dual basis for this field.
       
            .. warning::
                DO NOT USE, still under development.
        """

        # Compute a short list who's trace of their square is equal to 1
        first_round = []   
        for element in self.elements:
            if tr(element * element) == 1:
                first_round.append(element)
    
        for element in first_round:
            print(element)

        second_round = []

        # Of the remaining possible elements, compute traces and see 
        # if we can find n of them which equal 0
        for i in range(0, len(first_round)):
            traces = [tr(first_round[i] * first_round[j]) for j in range(0, len(first_round))] 
            if traces.count(0) == self.n:
                second_round.append(first_round[i])
                print(traces)

        print(second_round)

        return


    def to_poly(self):
        """ Switch back to representation in the polynomial basis. 
        """
        for el in self.elements:
            el.is_sdb = False
            el.sdb_coefs = []
        self.is_sdb = False


    def evaluate(self, coefs, argument):
        """ Evaluate a function, or curve on a finite field element.
            We consider here functions of the form
            
            .. math::
              
              f(\\alpha) = c_0 + c_1 \\alpha + \cdots + c_n \\alpha^n
            This function is primarily meant for use with the Curve class in
            my Balthasar package.
            Args:
                coefs (list): A set of coefficients for the curve, i.e.
                              :math:`[c_0, c_1, \ldots, c_n]`. These should
                              be a mix of integers and FieldElements.
                argument (FieldElement): The argument to the function, i.e.
                      the :math:`\\alpha` in :math:`f(\\alpha)`.
              
            Returns:
                The value of the function of the argument, taken over the
                finite field.
        """
        result = coefs[0] * self.elements[-1] 
        for coef_idx in range(1, len(coefs)):
            result += coefs[coef_idx] * pow(argument, coef_idx)
        return result


    def print(self):
        """ Print out all the useful information about a field."""
        
        print("--- Galois field information ---")
        print("p = " + str(self.p))

        if self.n > 1:
            print("n = " + str(self.n))

            print("\nIrreducible polynomial: ", end = "")
            if self.coefs[0] != 0:
                print(str(self.coefs[0]) + " + ", end ="")

            for i in range(1, len(self.coefs)):
                if self.coefs[i] != 0:
                    # Print coefficient if it's not 1
                    if self.coefs[i] != 1:
                        print(str(self.coefs[i]), end = "") 

                    # Print exponent value
                    print("x", end = "")
                    if i != 1:
                        print("^" + str(i), end = "")

                    if i != len(self.coefs) - 1: 
                        print(" + ", end = "")

        print("\nField elements:")
        for element in self.elements:
            element.print()


def tr(x):
    """ Wrapper trace function so the user can do tr(x) or x.trace()."""
    # Make sure x is a field element
    if type(x) is not FieldElement:
        print("Error, invalid argument to function 'tr'.")
        return None
    else:
        return x.tr()


def gchar(x):
    """ Wrapper so the user can do x.gchar() or gchar(x). """
    if type(x) is not FieldElement:
        print("Error, invalid argument to function 'gchar'.")
        return None
    else:
        return x.gchar()


def inv(x):
    """ Wrapper so the user can do x.inv() or inv(x) interchangeably."""
    # Make sure x is a field element
    if type(x) is not FieldElement:
        print("Error, invalid argument to function 'inv'.")
        return None
    else:
        return x.inv()



import math

class FieldElement():
    """ Class for an element in a finite field.
        Args:
            p (int): The prime order of the field this element is in.
            n (int): The degree of the field extension for the field this 
                     element is in.
            exp_coefs (list): The set of expansion coefficients of this element
                              in terms of some basis.
            field_list (list of FieldElements) - A copy of the list of all 
                     elements in the GaloisField this element is in. Empty when
                     the elements are initially constructed, and filled later
                     by the constructor of GaloisField. Inclusion 
                     of this parameter is not ideal, but greatly simplifies the 
                     implementation of many of the arithmetic operations 
                     (notably multiplication, inverse, exponentiation),
                     specifically for power of prime fields. With field_list, 
                     each element always knows its position, or power of the 
                     primitive element in the field. I'm not proud of it being
                     implemented this way, but this is the best I can do now.
        Attributes:
            p (int): The prime order of the field this element is in.
            n (int): The degree of the field extension for the field this 
                     element is in.
            dim (int): The dimension of the field, :math:`p^n`.
            prim_power (int): This element represented as a power of the 
                              primitive element of the field.
            exp_coefs (list): The set of expansion coefficients of this element
                              in terms of the polynomial basis.
            is_sdb (bool): Indicates whether this element is expressed in
                           the self-dual basis or not. False by default.
            sdb_coefs (list): The set of expansion coefficients in the self-
                              dual basis. Empty by default.
            str_rep (string): A representation of the exp_coefs as a string.
            field_list (list of FieldElements) - A copy of the list of all 
                     elements in the GaloisField this element is in.
    """

    def __init__(self, p, n, exp_coefs, field_list = [], is_sdb = False, sdb_field_list = [], sdb_coefs = []):
        self.p = p
        self.n = n
        self.dim = int(math.pow(p, n))

        # Set the expansion coefficients.
        # If we're in a prime field, the basis is 1, and
        # the coefficient is just the value
        self.exp_coefs = exp_coefs
        self.str_rep = ",".join([str(x) for x in exp_coefs])

        # These parameters initially gets set by the GaloisField constructor 
        # after ALL the field elements have been created. This is set only
        # for power of prime fields. However, when we perform operations on
        # elements such as addition, multiplication, we will need to 
        self.field_list = field_list

        if len(field_list) != 0:
            self.prim_power = self.field_list.index(self.str_rep)
        else:
            self.prim_power = -1

        # Prim power doesn't really make sense for prime, but set it here
        # so that we can make the rest of the code more general
        if self.n == 1:
            self.prim_power = exp_coefs[0]

        # These parameters will be something other than their default value
        # only if the to_sdb function is called on the GaloisField.
        self.is_sdb = is_sdb
        self.sdb_field_list = sdb_field_list
        self.sdb_coefs = sdb_coefs 

        # Reset the sdb coefficients after an operation if need be.
        if self.is_sdb:
            self.sdb_coefs = [int(x) for x in self.sdb_field_list[self.prim_power].split(',')] 
      

    def __add__(self, el):
        """ Addition.
            Args:
                el (FieldElement): A FieldElement to add to this one.
            Returns:
                A FieldElement which is this element + el. For prime fields
                this is simply addition modulo :math:`p`, for power-of-prime
                fields we must add using the exp_coefs.
        """
        # Make sure we're in the same field!
        if (self.p != el.p) or (self.n != el.n):
            print("Error, cannot add elements from different fields!")
            return None

        # Prime case
        if self.n == 1:
            return FieldElement(self.p, self.n, [(self.prim_power + el.prim_power) % self.p])
        else: # Power of prime case
            # Coefficients simply add modulo p
            new_coefs = [(self.exp_coefs[i] + el.exp_coefs[i]) % self.p for i in range(0, self.n)]
            return FieldElement(self.p, self.n, new_coefs, self.field_list, self.is_sdb, self.sdb_field_list)


    def __radd__(self, el):
        """ Add a field element to the left of this one. 
        
            Addition in finite fields is commutative so this works just like
            the normal add. This is implemented so we can use 'sum' 
            over lists of FieldElement.
        """
        return self + el
    

    def __sub__(self, el):
        """ Addition.
            Args:
                el (FieldElement): A FieldElement to subtract from this one.
            Returns:
                A FieldElement which is this element - el. For prime fields
                this is simply subtraction modulo :math:`p`, for power-of-prime
                fields we must subtract using the exp_coefs.
        """
        # Make sure we're in the same field!
        if (self.p != el.p) or (self.n != el.n):
            print("Error, cannot subtract elements from different fields!")
            return None

        # Prime case
        if self.n == 1:
            return FieldElement(self.p, self.n, [(self.prim_power - el.prim_power) % self.p])
        else:  # Power of prime case
            # Coefficients subtract modulo p
            new_coefs = [(self.exp_coefs[i] - el.exp_coefs[i]) % self.p for i in range(0, self.n)]
            return FieldElement(self.p, self.n, new_coefs, self.field_list, self.is_sdb, self.sdb_field_list)


    def __mul__(self, el):
        """ Multiplication.
            Args: 
                el (int or FieldElement): An element to multiply with this one.
                      Can also pass an integer value.
            Returns:
                This element * el. For prime fields, this amounts to simple
                multiplication modulo :math:`p`. For power of primes, this is
                where the ugly field_list comes in handy. We can compute the
                new power of the primitive element by adding together this one
                and the one from el; we then use field_list to find the 
                corresponding FieldElement and return it.
        """
        # Multiplication by a constant (must be on the right!)
        if isinstance(el, int):
            return FieldElement(self.p, self.n, [(el * exp_coef) % self.p for exp_coef in self.exp_coefs] , self.field_list, self.is_sdb, self.sdb_field_list)

        # Multiplication by another FieldElement
        elif isinstance(el, FieldElement):
            # Make sure we're in the same field!
            if (self.p != el.p) or (self.n != el.n):
                print("Error, cannot multiply elements from different fields!")
                return None

            # Prime case
            if self.n == 1:
                return FieldElement(self.p, self.n, [(self.prim_power * el.prim_power) % self.p])
            # Power of prime case
            else:
                # I stored the whole list of field elements in each element for a reason...
                # Now we can multiply really easily

                # Multiplying by 0, nothing to see here
                if el.prim_power == 0 or self.prim_power == 0: 
                    zeros = [0] * self.n
                    return FieldElement(self.p, self.n, zeros, self.field_list, self.is_sdb, self.sdb_field_list)
                else:
                    new_exp = self.prim_power + el.prim_power # New exponent
                    # If the exponent calculated is outside the range of primitive element
                    # powers of the field, we need to wrap it around using the fact that
                    # the last field element is 1.
                    if new_exp > self.dim - 1: 
                        new_exp = ((new_exp - 1) % (self.dim - 1)) + 1
                    new_exp_coefs = [int(x) for x in self.field_list[new_exp].split(',')] 
                    return FieldElement(self.p, self.n, new_exp_coefs, self.field_list, self.is_sdb, self.sdb_field_list)
        else:
            raise TypeError("Unsupported operator")


    def __rmul__(self, el): # Implementing rmul so we can multiply on the left by integers
        """ Multiplication from the left. """
        return self * el
 

    def __truediv__(self, el):
        """ Division.
            In a Galois Field division is just multiplying by the inverse. By
            definition of a finite field, every element has a multiplicative
            inverse, except for 0.
            Args:
                An element to divide this one by.
            Returns:
                This element / el. Returns None if el = 0. 
        """
        if isinstance(el, FieldElement):
            if (self.p != el.p) or (self.n != el.n):
                print("Error, cannot divide elements from different fields.")

            # Prime
            if self.n == 1:
                if self.prim_power == 0:
                    print("Cannot divide by 0.")
                    return
            # Power of prime
            else:
                if self.field_list.index(self.str_rep) == 0:
                    print("Cannot divide by 0.")
                    return
            # Actually do the division 
            return self * el.inv()


    # Operations with assignment
    def __iadd__(self, el):
        """ Addition with assignment. """
        return self + el


    def __isub__(self, el):
        """ Subtraction with assignment. """
        return self - el


    def __imul__(self, el):
        """ Multiplication with assignment. """
        return self * el


    def __itruediv__(self, el):
        """ Division with assignment. """
        return self / el


    def __pow__(self, exponent):
        """ Exponentiation.
            Args:
                exponent (int): Something to exponentiate this element by.
            Returns:
                This element to the power of exponent. Just the normal power
                modulo p for primes. For power-of-primes, we define that the
                power of any element to 0 is the 0 element, and *not* 1.
        """
        # Prime case
        if self.n == 1:
            return FieldElement(self.p, self.n, [int(math.pow(self.prim_power, exponent)) % self.p])
        # Power of prime case
        else:
            new_coefs = []
            # 0, and any element to the 0 is 0 by convention 
            if self.prim_power == 0 or exponent == 0: 
                new_coefs = [int(x) for x in self.field_list[0].split(',')]
            else:
                new_exp = self.prim_power * exponent
                if new_exp > self.dim - 1:
                    new_exp = ((new_exp - 1) % (self.dim - 1)) + 1
                new_coefs = [int(x) for x in self.field_list[new_exp].split(',')] 
            return FieldElement(self.p, self.n, new_coefs, self.field_list, self.is_sdb, self.sdb_field_list)
            

    def __eq__(self, el):
        """ Test equality of two field elements.
            
            Args:
                el (FieldElement): An element to compare with.
            Returns:
                True if the field dimensions (:math:`p`, :math:`n`) are the 
                same, the basis expansions are the same, and the list of 
                field elements is the same. False otherwise.
        """
        if (self.p != el.p) or (self.n != el.n):
            return False
        if self.exp_coefs != el.exp_coefs:
            return False
        if self.field_list != el.field_list:
            return False
        return True


    def __lt__(self, el):
        """ Implement a 'natural' ordering for field elements.
            For prime fields, this is simply the ordering of natural numbers.
            For power of primes, turn the coefficient lists into binary
            strings, and order them this way. Doing this to allow for
            Wigner functions to be plotted 'in order' in Balthasar.
            Args:
                el (FieldElement): An element to compare with.
            Returns:
                True if this element is 'less' by the conditions defined above.
                False otherwise.
        """
        if self.n == 1: 
            if self.prim_power < el.prim_power:
                return True
            else:
                return False
        else:
            # If there is a sdb defined, use that, otherwise use exp_coefs
            if self.is_sdb:
                this_exp_str = [str(x) for x in self.sdb_coefs]
                that_exp_str = [str(x) for x in el.sdb_coefs]
                if "".join(this_exp_str) < "".join(that_exp_str):
                    return True
                else:
                    return False
            else:
                this_exp_str = [str(x) for x in self.exp_coefs]
                that_exp_str = [str(x) for x in el.exp_coefs]
                if "".join(this_exp_str) < "".join(that_exp_str):
                    return True
                else:
                    return False


    def __repr__(self):
        """ Make the field element get printed in the command line."""
        if self.n == 1:
            return str(self.prim_power)
        else:
            if self.is_sdb:
                return str(self.sdb_coefs)
            else:
                return str(self.exp_coefs)


    def __hash__(self):
        """ Make hashable so we can use these guys as dictionary keys."""
        return hash(repr(self))


    def inv(self):
        """ Compute the multiplicative inverse of a field element.
            Returns:
                The FieldElement that is the inverse of this one. All
                elements have a multiplicative inverse except for 0;
                if 0 is passed, prints error message and returns None.
            Note: The trace of an element can be invoked in two ways. One can
            do el.inv() or inv(el).
        """
        if self.n == 1: # Prime case - brute force :(
            if self.prim_power == 0:
                print("Error, 0 has no multiplicative inverse.")
                return

            for i in range(0, self.p):
                if (self.prim_power * i) % self.p == 1:
                    return FieldElement(self.p, self.n, [i])
        else: # Power of prime case
            if self.prim_power == 0:
                print("Error, 0 has no multiplicative inverse.")
                return 
            # Last element is always 1 which is it's own inverse
            elif self.prim_power == self.dim - 1:
                return self 
            # All other elements, find exponent which sums to dim - 1
            else:
                new_coefs = [int(x) for x in self.field_list[self.dim - self.prim_power - 1].split(',')]
                return FieldElement(self.p, self.n, new_coefs, self.field_list, self.is_sdb, self.sdb_field_list)


    def tr(self):
        """ Compute the trace of a field element. 
            The trace of an element is defined as 
            .. math ::
              \\text{tr}(\\alpha) = \\alpha + \\alpha^p + \\alpha^{p^2} + \cdots \\alpha^{p^{n - 1}}
        
            The trace of any element should be an element of the base field 
            GF(:math:`p`) for the power of prime case.
            Returns:
                The trace of this element, as expressed above, as an integer.
            Note: The trace of an element can be invoked in two ways. One can
            do el.tr() or tr(el).
        """
        s = self

        if self.n == 1:
            return self.prim_power
        else:
            for i in range(1, self.n):
                s = s + pow(self, pow(self.p, i))

        return s.exp_coefs[0]


    def gchar(self):
        """ Compute the character of a FieldElement. 
            We define our group character as 
            .. math::
                \chi({\\alpha}) = \omega_{p}^{\\text{tr}({\\alpha})} 
            where :math:`\omega_p` is the :math:`p^\\text{th}` root of unity.
            Returns:
                :math:`\chi(\\alpha)` as defined above. For fields with 
                characteristic 2, this is an integer. For fields extended from
                odd primes, it is a pthRootOfUnity.
            Note: The trace of an element can be invoked in two ways. One can
            do el.gchar() or gchar(el).
        """
        if self.p == 2:
            return ((-1) **  self.tr())
        else:
            return pthRootOfUnity(self.p, self.tr())
            

    def print(self):
        """ Print out information about the element."""
        if self.n == 1:
            print(self.prim_power)
        else:
            if self.is_sdb:
                print(self.sdb_coefs)
            else:
                print(self.exp_coefs)



class pthRootOfUnity():
    """ Class to hold :math:`p^{\\text{th}}` roots of unity symbolically over finite fields. 
        
        These guys have the form
        .. math::
            \omega_p = \exp \left( \\frac{2 \\pi i}{p} \\right)
        where :math:`p` is the characteristic of the field.
        They can be evaluated both symbolically and also exactly by 
        explicitly computing the value above using numpy. 
        Here we'll implement only what we need: exponents and multiplication.
        Args:
            p (int): The characteristic of some finite field.
            e (int): An exponent, if desired (default is 1).
        An object of type pthRootOfUnity has the following attributes.
        Attributes:
            p (int): Which root of unity we are dealing with.
            e (int): An additional exponent representing a power, i.e. 
                     :math:`\omega_p ^ e`.
    """

    def __init__(self, p, e = 1):
        # TODO some sort of primality check
        if p < 2:
            print ("Error, please enter a prime to make roots of unity.")
        else:
            self.p = p
            self.e = e


    def __mul__(self, op):
        """ Multiply two roots of unity. Roots of unity are cyclic, i.e.
            :math:`\omega_p ^ p = 1`.
            Args:
                op (pthRootOfUnity): What to multiply by, say, 
                                     :math:`\omega_p^{e^\prime}`.
            Returns: 
                :math:`\omega_p^{e} \cdot \omega_p^{e^\prime}`.
        """
        if self.p != op.p:
            print("Error, cannot multiply roots of unity from different primes.")
            return 

        new_exp = (self.e + op.e) % self.p
        return pthRootOfUnity(self.p, new_exp)  


    def __imul__(self, op):
        """ Multiplication with assignment. """
        return self * op


    def __truediv__(self, op):
        """ Division. 
            Args:
                op (pthRootOfUnity): What to divide by, say, 
                                     :math:`\omega_p^{e^\prime}`.
            Returns: 
                :math:`\omega_p^{e} /  \omega_p^{e^\prime}`.
        """
        if self.p != op.p:
            print("Error, cannot divide roots of unity from different primes.")
            return 

        new_exp = (self.e - op.e) % self.p
        return pthRootOfUnity(self.p, new_exp)  


    def __itruediv__(self, op):
        """ Division with assignment. """
        return self / op


    def __pow__(self, ex):
        """ Exponentiation. 
        
            Args:
                ex (int): Some integer exponent.
            Returns:
                :math:`\omega_p^{e \cdot ex}`, where :math:`ex` is
                the exponent passed in as an argument.
        """
        new_exp = (self.e * ex) % self.p
        return pthRootOfUnity(self.p, new_exp)  

            
    def __eq__(self, op):
        """ Equality of two roots of unity.
            We consider two roots of unity equal if both their characteristics
            :math:`p` and exponents :math:`e` are equal.
        
            Args:
                op (pthRootOfUnity): A pth root to check equality with.
            Returns:
                True if the primes and exponents are the same; false otherwise.
                None if there is a type error.
        """
        if type(op) != pthRootOfUnity:
            print("Error, type cannot be compared with pthRootOfUnity.")
            return None

        if self.p != op.p:
            return False
        if self.e != op.e:
            return False
        return True


    def __repr__(self):
        """ Print pth roots of unity in the command line.
            We represent the :math:`\omega` using a simple 'w'.
        """ 
        return ("w^" + str(self.e))


    def eval(self):
        """ Use numpty to evaluate the actual number. Gross. 
        
            Returns: 
                The numerical value :math:`\exp \\left(\\frac{2 \\pi i \cdot e}{p} \\right)`.
        """
        return pow(np.exp(2 * np.pi * 1j / self.p), self.e)


    def print(self):
        """ Prints the pth root of unity. 
        """
        print("w^" + str(self.e))



def generate_field(p, primitive):
  primitive = primitive[::-1]
  print(primitive)

  ges_wrong_order = GaloisField(2, p, primitive).elements

  ges=[]
  for ge in ges_wrong_order: 
    ges.append(ge.exp_coefs[::-1])

  ges.insert(1, ges.pop())

  return ges

