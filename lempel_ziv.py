import math



def binary_repesenation(number, bits):
	"""
	Return the "number" in a binary string of "bits" length 
	"""
	rep = ""

	while bits > 0:
		bit_value = 2**(bits-1)
		if number/bit_value>=1:
			rep += "1"
			number -= bit_value
		else:
			rep += "0"
		bits -= 1

	return rep

def number_from_binary(binary_string):
	number = 0
	for i in range(len(binary_string)):
		s = binary_string[-(i+1)]

		number += int(s) * 2**i
	return number


def lempel_ziv_encode(binary_string):

	binary_string=str(binary_string)

	lamb = ""
	phrases = [lamb]

	# 1. Generate the phrases
	current_phrase = ""
	output_string = ""
	for s in binary_string:
		
		current_phrase += s
		if (current_phrase not in phrases):
			

			pointer_number = phrases.index(current_phrase[:-1])
			
			phrases_length = len(phrases)

			# Build a pointer
			if len(phrases) == 1:
				# Redundant pointer, so add no pointer
				pointer = ""

			else:
				# Pointers are only as long as they need to be, go up in bit length once we have too many phrases to count


				# this gives the right number of bits to add
				pointer_bits = next_pointer_bit_length(phrases_length)

				pointer = binary_repesenation(pointer_number, pointer_bits)


			phrases.append(current_phrase)

			output_string += pointer + current_phrase[-1]

			# Reset the current phrase
			current_phrase = ""

	if len(current_phrase) > 0:
		output_string += current_phrase

	return output_string, phrases


def next_pointer_bit_length(phrases_length):
	if phrases_length == 1:
		return 0
	return math.floor(math.log2(phrases_length-1)+1)



def lempel_ziv_decode(binary_string):

	lamb = ""
	phrases = [lamb]

	decoded_string = ""

	binary_string_length = len(binary_string)
	next_block_start_position = 0
	for s in binary_string:

		phrases_length = len(phrases)

		next_pointer_length = next_pointer_bit_length(phrases_length)
		next_block_length = next_pointer_length + 1

		# At the end of the string, jsut add the remaining symbols to the output
		if next_block_start_position + next_block_length > binary_string_length:
			decoded_string += binary_string[next_block_start_position:]
			break
		# If its the first block, there is a redundant pointer so handle it as a special case
		elif len(phrases) == 1:

			current_phrase = binary_string[0]
			decoded_string += current_phrase
			phrases.append(current_phrase) 
			next_block_start_position = 1
		else:

			current_block = binary_string[next_block_start_position:next_block_start_position+next_block_length]
			
			current_pointer = current_block[:-1]


			phrase_number = number_from_binary(current_pointer)
			


			current_phrase = phrases[phrase_number] + current_block[-1]
			phrases.append(current_phrase)

			decoded_string += current_phrase

			next_block_start_position += next_block_length


	return decoded_string






def test_lempel_ziv_encoding():

	start_string = "000000000000100000000000"
	encoded_string = lempel_ziv_encode(start_string)
	target_string = "010100110010110001100"
	if encoded_string == target_string:
		print("Test 1 passed")
	else:
		print("TEST FAILED!")

	reversed_string = lempel_ziv_decode(target_string)
	if reversed_string == start_string:
		print("Test 1 decoding passed")
	else:
		print("TEST FAILED!")


	start_string = "1011010100010"
	encoded_string = lempel_ziv_encode(start_string)
	target_string = "100011101100001000010" 
	if encoded_string == target_string:
		print("Test 2 passed")
	else:
		print("TEST FAILED!")

	reversed_string = lempel_ziv_decode(target_string)
	if reversed_string == start_string:
		print("Test 2 decoding passed")
	else:
		print("TEST FAILED!")

	start_string = "10110101000100"
	encoded_string = lempel_ziv_encode(start_string)
	target_string = "1000111011000010000100" 
	if encoded_string == target_string:
		print("Test 3 passed")
	else:
		print(encoded_string)
		print("TEST FAILED!")

	reversed_string = lempel_ziv_decode(target_string)
	if reversed_string == start_string:
		print("Test 3 decoding passed")
	else:
		print("TEST FAILED!")

	start_string = "Hi there john boy whats up"
	encoded_string = lempel_ziv_encode(start_string)
	print(encoded_string)

if __name__=="__main__":
	test_lempel_ziv_encoding()