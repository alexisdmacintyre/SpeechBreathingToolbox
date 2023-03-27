# SpeechBreathingToolbox
Tools for the automatic detection of speech-related inhalation events and characterisation of the speech respiratory cycle.

Please cite as:

MacIntyre, A. D. and Werner, R. (2023). An Automatic Method of Speech
Breathing Annotation. <i>Proceedings of the 34th Conference on Electronic
Speech Signal Processing (ESSV)</i>, Munich, DE.

**About**

Use `findBreaths()` to automatically detect speech breathing onsets and ends, then optionally corroborate the results to the acoustic speech envelope using `breathSpeechCompare()`. Finally, you can use `plotBreaths()` to visualise the breath onsets and ends in relation to the breath belt signal and, if you have it, the corresponding acoustic data. This plot can help you decide if you need to adjust the default paramaters (e.g., minimum inter-breath interval).

**Usage**

Run `example.m` with included example breath belta data and acoustic recording to see how the scripts work.

<img width="900" alt="Acoustic landmarks for speech rhythm analysis" src="https://user-images.githubusercontent.com/55560694/215552770-4264208e-aaf1-4a16-9365-4db1460a0b8a.png">

**Future Versions**

Plans include a GUI to adjust or reject inhalation events, as well as extract basic breath- and/or speech-related descriptive statistics. As the 'SpeechBreathingToolbox' is under development, reporting issues is appreciated.
