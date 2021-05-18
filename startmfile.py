import matlab.engine

print("Loading matlab engine...")
eng = matlab.engine.start_matlab()
eng.main(nargout=0)
