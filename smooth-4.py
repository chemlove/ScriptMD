# smooths a list of values using weighted moving average, up to N values
# in each direction
def smooth(data, N):
    datalen = len(data)
    data_smoothed = []
    for t in range(datalen):

	#T stands for width of window
	if t < N:
		T = min(t, datalen - t - 1) + 1
	else:
		T = N

        numer = T * data[t]
        denom = T

        for n in range(1,T):
            if t-n >= 0:
                numer += (T-n) * data[t-n]
                denom += (T-n)
            if t+n < datalen:
                numer += (T-n) * data[t+n]
                denom += (T-n)
        data_smoothed.append(numer / float(denom))
    return data_smoothed


