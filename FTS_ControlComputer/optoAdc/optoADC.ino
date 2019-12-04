/*
  Analog input, analog output, serial output
  Reads an analog input pin, maps the result to a range from 0 to 255 and uses
  the result to set the pulse width modulation (PWM) of an output pin.
  Also prints the results to the Serial Monitor.

*/

// These constants won't change. They're used to give names to the pins used:
const int analogInPin = A0;  // Analog input pin that the potentiometer is attached to
const int analogOutPin = 9; // Analog output pin (digital)

int sensorValue = 0;        // voltage signal from the opto-interrupter 
                            // Arduino references from 5V as a default, A value of 1024 corresponds to 5V signal


int outputValue = 0;        // value output to the PWM (analog out)
int sensorMin = 1023;
int sensorMax = 0;
unsigned long time;


void setup(){

  //calibrate first 5 seconds
  while(millis() < 5000){
    sensorValue = analogRead(analogInPin);
    if(sensorValue > sensorMax){
      sensorMax = sensorValue;
    }
    if(sensorValue < sensorMin){
      sensorMin = sensorValue;
    }
  }
  // initialize serial communications at 9600 bps:
  //change this to change the output rate
    Serial.begin(115200);
}
void loop() {
  // read the analog in value:
  sensorValue = analogRead(analogInPin);
  // map it to the range of the analog out:
  
  
  outputValue = map(sensorValue, sensorMin, sensorMax, 0, 255);
  // change the analog out value:
  analogWrite(analogOutPin, outputValue);

  // print the results to the Serial Monitor:
  time = millis();
 // Serial.print("sensor = ");
  Serial.print(sensorValue);
  Serial.print(" ");
  Serial.print(outputValue);
  Serial.print(" ");
  Serial.println(time);
 // Serial.println(" ");


  // wait 2 milliseconds before the next loop for the analog-to-digital
  // converter to settle after the last reading:
  //delay(2);
}
