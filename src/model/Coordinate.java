package model;

public class Coordinate {
	public double x;
	public double y;
	
	public Coordinate(double x, double y) {
		this.x = x;
		this.y = y;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("[" + x + ", " + y + "]");
		return sb.toString();
	}
}
